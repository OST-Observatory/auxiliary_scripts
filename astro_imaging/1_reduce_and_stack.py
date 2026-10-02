#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Astro imaging, step 1 of 2: reduce, select the best frames, register
    every filter onto ONE pixel grid, and write a weighted linear stack per
    filter.

    The result (``output/combined_filter_<F>.fit``, 32-bit float, e-/s,
    sigma-clipped weighted average) is the linear master for any later
    post-processing in Siril, PixInsight, GIMP, ... No stretching and no
    colour combination happens here, on purpose.

    What the script does with the parameters below:
      * bias / dark / flat calibration and optional cosmic-ray removal
      * FWHM, roundness, star count and sky background of every light frame
        (``output/frame_quality.ecsv`` and the FITS headers)
      * rejection of poor frames (threshold, best X %, roundness, ...);
        rejected frames are MOVED to ``output/rejected_lights/`` and kept
      * registration of all filters onto the grid of the globally sharpest
        frame (``shift_all=True``), so colour channels line up pixel-exactly
      * weighted stacking per filter, weights from FWHM / star count / noise
      * WCS solution on the stacks (for annotation and later plate matching)
      * QC plots: ``output/diagnostics/frame_quality/frame_quality_<F>.pdf``

    Iterate on the selection WITHOUT re-reducing: keep
    ``keep_aligned_lights = True`` and run ``2_restack.py``.

    All files can be given in one directory called 'raw_files'. Alternatively
    images can be sorted into the following directory structure:
        * Images of the object
        * Dark frames
        * Flatfields
        * Bias (optional)
    If they are sorted into directories the FITS header keywords will be
    checked for consistency.

    Images in sub folders will be recognized, but only one level is considered.

    Needs ost_photometry >= 0.6 (frame quality and weighted stacking).
"""

############################################################################
#                          Simple folder structure                         #
############################################################################
#   PLEASE NOTE: Either specify a single directory with all raw files at
#                this point, or use the directory structure shown below.
raw_files: str = "?"

############################################################################
#                              Individual folders                          #
############################################################################
#   Path to the bias -- If set to '?', bias exposures are not used.
bias: str = "?"

#   Path to the darks
darks: str = "?"

#   Path to the flats
flats: str = "?"

#   Path to the images
images: str = "?"


############################################################################
#                Additional options: only edit if necessary                #
############################################################################

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir: str = "output/"

#   Only reduce images whose FITS ``OBJECT`` keyword equals this name.
#   ``None`` uses every light frame.
target_name: str | None = None

#   Remove cosmic rays (L.A.Cosmic): True, False or "auto". "auto" leaves
#   them to the sigma-clipped stack for filters with at least
#   ``cosmic_ray_auto_min_frames`` frames and removes them otherwise. Avoid
#   True for strongly undersampled stars (FWHM < ~2 px), where L.A.Cosmic
#   can nibble at star cores.
rm_cosmic_rays: bool | str = "auto"
# rm_cosmic_rays: bool | str = True
cosmic_ray_auto_min_frames: int = 7

#   Gain and read noise: "catalog" (header / camera catalog, read noise
#   scaled to the binned pixel) or "measured" (from the bias and flat pairs
#   of this data set; falls back to the catalog if that is not possible).
camera_noise_source: str = "catalog"
# camera_noise_source: str = "measured"

#   How the camera bins, for the catalog read noise: "auto" (camera
#   catalog), "digital" (CMOS: read noise x sqrt(xbin * ybin)) or "charge"
#   (CCD on-chip binning).
binning_mode: str = "auto"

#   Tolerance between science and dark exposure times in s
exposure_time_tolerance: float = 5.0

#   Tolerance between the camera chip temperatures of the images
temperature_tolerance: float = 5.0

#   Number of cores used for multiprocessing (``None`` = half the CPUs)
n_cores_multiprocessing: int | None = None

#   Alignment backend. ``aa_true`` (astroalign similarity: shift, rotation,
#   scale; sub-pixel) is the imaging default. ``wcs`` reprojects onto the
#   sky grid of the reference and is the fallback for large dithers, sparse
#   fields or when astroalign fails (needs ASTAP on PATH).
shift_method: str = "aa_true"
# shift_method: str = "wcs"

############################################################################
#                       Frame quality and selection                        #
############################################################################
#   Every criterion is optional (``None`` = off). All set criteria must pass
#   (logical AND), evaluated per filter. Frames that show no stars at all
#   (clouds) are rejected as soon as any criterion is active.

#   Reject frames whose FWHM exceeds this value in pixels ...
fwhm_max_px: float | None = None
#   ... or in arcsec (needs FOCALLEN and XPIXSZ in the header, or a WCS).
#   Set only one of the two.
fwhm_max_arcsec: float | None = None

#   Keep only the sharpest fraction of the frames that pass the thresholds
#   (Siril "best X %"): 0.8 keeps the best 80 %. ``None`` keeps all.
best_fraction: float | None = 0.8

#   Rank for ``best_fraction``: ``fwhm_px`` (sharpness only) or
#   ``fwhm_weighted`` (Siril wFWHM: FWHM penalised by a low star count,
#   i.e. thin clouds count against a frame).
rank_by: str = "fwhm_px"

#   Reject frames whose FWHM lies more than k robust sigmas (MAD) above the
#   median of the filter, e.g. 3.0. ``None`` = off.
fwhm_sigma_clip: float | None = None

#   Maximum IRAF roundness (0 = round, 1 = a line). Tracking errors, wind
#   gusts and focus drift show up here. 0.3 is a good imaging default.
roundness_max: float | None = 0.3

#   Minimum number of detected stars (transparency, thin clouds).
n_stars_min: int | None = None

#   Maximum sky background in e-/s/pixel (moon, twilight, haze).
background_max: float | None = None

#   Never keep fewer frames than this per filter; the best rejected frames
#   are restored if the selection would cut deeper.
min_frames: int = 5

############################################################################
#                                Stacking                                  #
############################################################################
#   Combine method: ``average`` (sigma-clipped, supports weights), ``median``
#   (ignores weights), ``sum``.
stack_method: str = "average"

#   Per-frame weights for the stack: ``none``, ``fwhm`` ((median/FWHM)^2,
#   sharp frames dominate), ``n_stars`` (transparency), ``noise``
#   ((median RMS / RMS)^2). Weights are normalised per filter and clipped to
#   0.1-10 so a single frame cannot dominate.
stack_weighting: str = "fwhm"

#   Alignment reference: ``best_fwhm`` uses the sharpest kept frame of the
#   whole data set as the common grid; ``first`` uses the first frame.
reference_image_selection: str = "best_fwhm"

#   Or force a specific frame (index after sorting by time, over all
#   filters). An explicit index always wins over the selection above.
reference_image_index: int | None = None

#   Keep ``output/aligned_lights/`` (registered single frames). Needed for
#   ``2_restack.py`` and useful for drizzle in Siril. Costs disk space.
keep_aligned_lights: bool = True

#   Solve a WCS on the stacked images (ASTAP by default). Useful for later
#   annotation / plate matching; disable if ASTAP is not installed.
find_wcs: bool = True
wcs_method: str = "astap"

#   Data type for the combine step. ``None`` = float64 (safest);
#   ``"float32"`` halves the memory for very large sensors.
dtype: str | None = None

#   Pixels to trim from the frame edges (e.g. overscan).
trim_x_start: int = 0
trim_x_end: int = 0
trim_y_start: int = 0
trim_y_end: int = 0

############################################################################
#                               Libraries                                  #
############################################################################

import tempfile
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from astropy import log

log.setLevel("ERROR")

from ost_photometry import style
from ost_photometry.reduce import redu, utilities


############################################################################
#                         Helper: build the selection                      #
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


############################################################################
#                                  Main                                    #
############################################################################

if __name__ == "__main__":
    #   Set start time
    start_time = time.time()

    ###
    #   Prepare directories and make checks
    #
    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    #   Prepare directories
    raw_files = utilities.prepare_reduction(
        output_dir,
        bias,
        darks,
        flats,
        images,
        raw_files,
        temp_dir,
    )

    ###
    #   Reduce, select, register onto one grid, stack per filter
    #
    redu.reduce_main(
        raw_files,
        output_dir,
        rm_cosmic_rays=rm_cosmic_rays,
        cosmic_ray_auto_min_frames=cosmic_ray_auto_min_frames,
        camera_noise_source=camera_noise_source,
        binning_mode=binning_mode,
        exposure_time_tolerance=exposure_time_tolerance,
        temperature_tolerance=temperature_tolerance,
        target_name=target_name,
        n_cores_multiprocessing=n_cores_multiprocessing,
        shift_method=shift_method,
        #   Imaging: every filter on the same pixel grid, stack per filter
        shift_all=True,
        stack_images=True,
        stack_method=stack_method,
        #   Frame quality, selection, weights, reference
        measure_frame_quality=True,
        frame_selection=build_frame_selection(),
        stack_weighting=stack_weighting,
        reference_image_selection=reference_image_selection,
        reference_image_index=reference_image_index,
        keep_aligned_lights=keep_aligned_lights,
        #   WCS on the stacks only
        find_wcs=find_wcs,
        wcs_method=wcs_method,
        find_wcs_of_all_images=False,
        dtype=dtype,
        trim_x_start=trim_x_start,
        trim_x_end=trim_x_end,
        trim_y_start=trim_y_start,
        trim_y_end=trim_y_end,
    )

    out = Path(output_dir)
    stacks = sorted(p.name for p in out.glob("combined_filter_*.fit"))
    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print(f"   Stacks:          {', '.join(stacks) if stacks else 'none written'}")
    print(f"   Frame quality:   {out / 'frame_quality.ecsv'}")
    print(f"   Rejected frames: {out / 'rejected_lights'}")
    if keep_aligned_lights:
        print(f"   Aligned frames:  {out / 'aligned_lights'}  (input for 2_restack.py)")
    print(f"   QC plots:        {out / 'diagnostics' / 'frame_quality'}")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
