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

n_cores_multiprocessing: int | None = None

############################################################################
#                            Selection / stacking                          #
############################################################################
#   Same meaning as in astro_imaging/1_reduce_and_stack.py; ``None`` = off.
fwhm_max_px: float | None = None
best_fraction: float | None = None
roundness_max: float | None = None
n_stars_min: int | None = None
min_frames: int = 1

stack_method: str = "average"
#   none | fwhm | n_stars | noise
stack_weighting: str = "fwhm"

#   "wcs" (rotation / different cameras; uses ASTAP when needed) or "aa_true".
shift_method: str = "wcs"

#   "separate": one stack per camera (photometry); "combine": additionally a
#   noise-weighted stack over all cameras per filter (imaging).
camera_combination: str = "separate"

#   Keep the registered single frames (``stacks/<target>/aligned_lights/``).
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


def frame_selection() -> dict:
    selection: dict = {"min_frames": min_frames}
    for key, value in (("fwhm_max", fwhm_max_px), ("best_fraction", best_fraction),
                       ("roundness_max", roundness_max), ("n_stars_min", n_stars_min)):
        if value is not None:
            selection[key] = value
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
            frame_selection=frame_selection(),
            stack_weighting=stack_weighting,
            stack_method=stack_method,
            shift_method=shift_method,
            camera_combination=camera_combination,
            keep_aligned_lights=keep_aligned_lights,
            keep_reduced_lights=keep_reduced_lights,
            min_frames=min_frames,
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
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
