#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Archive pipeline, optional step 4: re-stack the registered frames of
    every target with a different frame selection and/or weighting, WITHOUT
    reducing and registering again.

    Reads ``<output_dir>/stacks/<target>/frame_quality.ecsv`` and the
    registered frames in ``<output_dir>/stacks/<target>/aligned_lights/``
    written by ``3_reduce_and_stack.py`` (with ``keep_aligned_lights =
    True``), applies the selection below per target, camera and filter, and
    writes new stacks, quality tables and QC plots to ``restack_dir``.
    Frames are only excluded here, never moved.

    Loosening the selection can only bring back frames that step 3
    registered: frames rejected there are available only if step 3 ran with
    ``align_rejected = True``.

    Typical loop: look at ``stacks/<target>/diagnostics/frame_quality/*.pdf``,
    tighten or loosen ``best_fraction`` / ``fwhm_max_*`` / ``roundness_max``,
    run this script, compare the stacks.

    Output in ``restack_dir``:
      * ``<target>/<camera>/combined_filter_<F>.fit``
      * ``<target>/combined_filter_<F>.fit``  (camera_combination = "combine")
      * ``<target>/frame_quality.ecsv`` and QC plots
      * ``summary.ecsv``

    Needs ost_photometry >= 0.6.
"""

############################################################################
#                                  Input                                   #
############################################################################
#   Working directory of steps 1 and 2 (for the target list of the plan).
workspace: str = "archive_workspace/"

#   Output directory of ``3_reduce_and_stack.py``.
output_dir: str = "output/"

#   Where the new stacks go (created if missing).
restack_dir: str = "output/restack/"

#   Re-stack only these targets (names or ids from the plan); ``None`` = all
#   targets with ``stack: true``.
targets: list[str] | None = None

############################################################################
#                       Frame quality and selection                        #
############################################################################
#   Same meaning as in ``3_reduce_and_stack.py``; ``None`` = off. Evaluated
#   per target, camera and filter.
fwhm_max_px: float | None = None
fwhm_max_arcsec: float | None = None
best_fraction: float | None = 0.8
rank_by: str = "fwhm_px"
fwhm_sigma_clip: float | None = None
roundness_max: float | None = 0.3
n_stars_min: int | None = None
background_max: float | None = None
min_frames: int = 5

############################################################################
#                                Stacking                                  #
############################################################################
#   ``average`` (sigma-clipped, supports weights), ``median``, ``sum``.
stack_method: str = "average"

#   ``none``, ``fwhm``, ``n_stars`` or ``noise``.
stack_weighting: str = "fwhm"

#   Groups with fewer kept frames are not stacked.
min_frames_per_stack: int = 1

#   "separate": one stack per camera; "combine": additionally a
#   noise-weighted stack over all cameras per filter.
camera_combination: str = "separate"

#   Data type for the combine step. ``None`` = float64.
dtype: str | None = None

############################################################################
#                               Libraries                                  #
############################################################################

import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from astropy import log

log.setLevel("ERROR")

from ost_photometry import style
from ost_photometry.reduce.grouping.plan import load_plan
from ost_photometry.reduce.workflow.combine import StackSettings, restack_planned

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
    target_dir = Path(restack_dir).expanduser()
    if not (out / "stacks").is_dir():
        raise SystemExit(
            f"{style.Bcolors.FAIL}No {out / 'stacks'}. Run 3_reduce_and_stack.py first "
            f"(with keep_aligned_lights = True).{style.Bcolors.ENDC}"
        )
    plan = load_plan(work / "calibration_plan.yaml")

    summary = restack_planned(
        plan,
        out,
        target_dir,
        StackSettings(
            frame_selection=build_frame_selection(),
            stack_weighting=stack_weighting,
            stack_method=stack_method,
            camera_combination=camera_combination,
            min_frames=min_frames_per_stack,
            dtype=dtype,
        ),
        targets=targets,
        log=lambda text: print("   " + text, flush=True),
    )
    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    for row in summary:
        print(f"   {row['target_name']:20s} {row['camera']:12s} {row['filter']:8s} "
              f"{row['n_images']:4d} frames {row['exposure_s'] / 60:7.1f} min  {row['path']}")
    print(f"   Summary: {target_dir / 'summary.ecsv'}")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
