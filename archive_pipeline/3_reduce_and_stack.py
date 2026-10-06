#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Archive pipeline, step 3 of 3: reduce the calibration plan and stack
    every target.

      * one bias / dark / flat master per calibration group,
      * lights reduced per unit (mount session x electronic setup) with
        exactly the masters the plan assigns,
      * per target: frame quality, selection, registration of all filters
        and cameras onto one grid (sharpest frame of the finest pixel
        scale), weighted stack per camera and filter over all nights,
      * optionally a noise-weighted combination of the cameras.

    Re-stack with another selection / weighting without reducing and
    registering again: ``4_restack.py`` (needs ``keep_aligned_lights``).

    Output in ``output_dir``:
      * ``masters/<master_id>/``             bias / dark / flat masters
      * ``reduced/<unit_id>/``               reduced lights (e-/s)
      * ``stacks/<target>/<camera>/combined_filter_<F>.fit``
      * ``stacks/<target>/combined_filter_<F>.fit``  (camera_combination = "combine")
      * ``stacks/<target>/frame_quality.ecsv`` and QC plots
      * ``stacks/summary.ecsv``, ``reduction_report.ecsv``

    Needs ost_photometry >= 0.6 (and ASTAP for shift_method = "wcs" when the
    frames carry no WCS yet).
"""

############################################################################
#                                  Input                                   #
############################################################################
#   Working directory of steps 1 and 2.
workspace: str = "archive_workspace/"

#   Where masters, reduced frames and stacks go.
output_dir: str = "output/"

#   Stack only these targets (names or ids from the plan); ``None`` = all
#   targets with ``stack: true``.
targets: list[str] | None = None

############################################################################
#                                Reduction                                 #
############################################################################
#   Cosmic-ray removal (L.A.Cosmic): True, False or "auto". "auto" removes
#   cosmics only from frames whose stack (target x camera x filter) has fewer
#   than ``cosmic_ray_auto_min_frames`` frames; larger stacks reject them by
#   sigma clipping, without eating into the sky noise.
rm_cosmic_rays: bool | str = "auto"
cosmic_ray_auto_min_frames: int = 7

#   Gain and read noise: "catalog" (header / camera catalog, read noise
#   scaled to the binned pixel) or "measured" (per electronic setup from bias
#   and flat pairs; falls back to the catalog where that is not possible).
camera_noise_source: str = "catalog"
# camera_noise_source: str = "measured"

#   How the camera bins, for the catalog read noise: "auto" (camera
#   catalog), "digital" (CMOS: read noise x sqrt(xbin * ybin)) or "charge"
#   (CCD on-chip binning).
binning_mode: str = "auto"

#   Electronics overrides (``None`` = see camera_noise_source). ``read_noise``
#   is per pixel of the images, i.e. per binned pixel.
gain: float | None = None
read_noise: float | None = None

#   Reuse masters of a previous run.
reuse_masters: bool = True

#   Lights the plan marks 'incomplete' (darks or flats missing) are skipped.
#   True reduces them anyway; to release single units use
#   overrides.force_units in calibration_plan.yaml and rerun step 2.
reduce_incomplete: bool = False

n_cores_multiprocessing: int | None = None

############################################################################
#                       Frame quality and selection                        #
############################################################################
#   Every criterion is optional (``None`` = off). All set criteria must pass
#   (logical AND), evaluated per target, camera and filter. Frames that show
#   no stars at all (clouds) are rejected as soon as any criterion is active.

#   Reject frames whose FWHM exceeds this value in pixels ...
fwhm_max_px: float | None = None
#   ... or in arcsec (needs FOCALLEN and XPIXSZ in the header, or a WCS;
#   the better choice when several cameras / binnings are combined).
#   Set only one of the two.
fwhm_max_arcsec: float | None = None

#   Keep only the sharpest fraction of the frames that pass the thresholds
#   (Siril "best X %"): 0.8 keeps the best 80 %. ``None`` keeps all.
best_fraction: float | None = None

#   Rank for ``best_fraction``: ``fwhm_px`` (sharpness only) or
#   ``fwhm_weighted`` (Siril wFWHM: FWHM penalised by a low star count,
#   i.e. thin clouds count against a frame).
rank_by: str = "fwhm_px"

#   Reject frames whose FWHM lies more than k robust sigmas (MAD) above the
#   median of the group, e.g. 3.0. ``None`` = off.
fwhm_sigma_clip: float | None = None

#   Maximum IRAF roundness (0 = round, 1 = a line). Tracking errors, wind
#   gusts and focus drift show up here. 0.3 is a good imaging default.
roundness_max: float | None = None

#   Minimum number of detected stars (transparency, thin clouds).
n_stars_min: int | None = None

#   Maximum sky background in e-/s/pixel (moon, twilight, haze).
background_max: float | None = None

#   Never keep fewer frames than this per group; the best rejected frames
#   are restored if the selection would cut deeper.
min_frames: int = 5

#   Register also the rejected frames (they are not stacked) so that
#   4_restack.py can loosen the selection later. Costs time and disk space.
align_rejected: bool = False

############################################################################
#                                Stacking                                  #
############################################################################
#   Combine method: ``average`` (sigma-clipped, supports weights), ``median``
#   (ignores weights), ``sum``.
stack_method: str = "average"

#   Per-frame weights for the stack: ``none``, ``fwhm`` ((median/FWHM)^2,
#   sharp frames dominate; arcsec when pixel scales are mixed), ``n_stars``
#   (transparency), ``noise`` ((median RMS / RMS)^2). Weights are normalised
#   per group and clipped to 0.1-10 so a single frame cannot dominate.
stack_weighting: str = "fwhm"

#   Groups (target x camera x filter) with fewer kept frames are not stacked.
min_frames_per_stack: int = 1

#   "wcs" (rotation / different cameras; uses ASTAP when needed) or "aa_true".
shift_method: str = "wcs"

#   "separate": one stack per camera (photometry); "combine": additionally a
#   noise-weighted stack over all cameras per filter (imaging).
camera_combination: str = "separate"

#   Keep the registered single frames (``stacks/<target>/aligned_lights/``);
#   needed for 4_restack.py.
#   With True, the reduced frames (``reduced/<unit>/``) of every aligned frame
#   are deleted after stacking unless ``keep_reduced_lights`` is True, so each
#   frame is stored once.
keep_aligned_lights: bool = True
keep_reduced_lights: bool = False

#   Floating type of masters, reduced / aligned frames and stacks: "float32"
#   (half the disk space, rounding far below the pixel noise) or "float64".
storage_dtype: str = "float32"

############################################################################
#                               Libraries                                  #
############################################################################

import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from astropy import log

log.setLevel("ERROR")

from astropy.table import Table

from ost_photometry import style
from ost_photometry.reduce.grouping.plan import load_plan
from ost_photometry.reduce.workflow.combine import StackSettings, stack_planned
from ost_photometry.reduce.workflow.groups import ReductionSettings, reduce_planned

############################################################################
#                                  Main                                    #
############################################################################


def build_frame_selection() -> dict:
    """Translate the parameter block into a ``frame_selection`` mapping."""
    if fwhm_max_px is not None and fwhm_max_arcsec is not None:
        raise ValueError("Set either fwhm_max_px or fwhm_max_arcsec, not both.")
    selection: dict = {"min_frames": min_frames, "rank_by": rank_by}
    if fwhm_max_px is not None:
        selection["fwhm_max"] = float(fwhm_max_px)
        selection["fwhm_unit"] = "px"
    elif fwhm_max_arcsec is not None:
        selection["fwhm_max"] = float(fwhm_max_arcsec)
        selection["fwhm_unit"] = "arcsec"
    if best_fraction is not None:
        selection["best_fraction"] = float(best_fraction)
    if fwhm_sigma_clip is not None:
        selection["fwhm_sigma_clip"] = float(fwhm_sigma_clip)
    if roundness_max is not None:
        selection["roundness_max"] = float(roundness_max)
    if n_stars_min is not None:
        selection["n_stars_min"] = int(n_stars_min)
    if background_max is not None:
        selection["background_max"] = float(background_max)
    return selection


if __name__ == "__main__":
    start_time = time.time()
    work = Path(workspace).expanduser()
    out = Path(output_dir).expanduser()
    plan = load_plan(work / "calibration_plan.yaml")
    frames = Table.read(work / "calibration_groups.ecsv", format="ascii.ecsv")

    reports, reduction = reduce_planned(
        plan,
        frames,
        out,
        ReductionSettings(
            rm_cosmic_rays=rm_cosmic_rays,
            cosmic_ray_auto_min_frames=cosmic_ray_auto_min_frames,
            camera_noise_source=camera_noise_source,
            binning_mode=binning_mode,
            gain=gain,
            read_noise=read_noise,
            n_cores_multiprocessing=n_cores_multiprocessing,
            reuse_masters=reuse_masters,
            storage_dtype=storage_dtype,
            reduce_incomplete=reduce_incomplete,
        ),
        log=lambda text: print("   " + text, flush=True),
    )
    for unit_id, report in reports.items():
        for note in dict.fromkeys(report.notes):
            print(f"   {unit_id}: {note}")

    summary = stack_planned(
        plan,
        frames,
        reduction,
        out,
        StackSettings(
            frame_selection=build_frame_selection(),
            stack_weighting=stack_weighting,
            stack_method=stack_method,
            shift_method=shift_method,
            camera_combination=camera_combination,
            keep_aligned_lights=keep_aligned_lights,
            keep_reduced_lights=keep_reduced_lights,
            align_rejected=align_rejected,
            min_frames=min_frames_per_stack,
            n_cores_multiprocessing=n_cores_multiprocessing,
        ),
        targets=targets,
        log=lambda text: print("   " + text, flush=True),
    )
    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    for row in summary:
        print(f"   {row['target_name']:20s} {row['camera']:12s} {row['filter']:8s} "
              f"{row['n_images']:4d} frames {row['exposure_s'] / 60:7.1f} min  {row['path']}")
    print(f"   Summary: {out / 'stacks' / 'summary.ecsv'}")
    if keep_aligned_lights:
        print("   Next (optional): 4_restack.py with another selection / weighting")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
