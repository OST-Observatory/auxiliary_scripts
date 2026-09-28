#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Astro imaging, step 2 of 2: re-stack the registered frames with a
    different frame selection and/or weighting, WITHOUT running the
    reduction and registration again.

    Reads ``<output_dir>/frame_quality.ecsv`` and the registered frames in
    ``<output_dir>/aligned_lights/`` written by ``1_reduce_and_stack.py``
    (with ``keep_aligned_lights = True``), applies the selection below and
    writes new stacks, a new quality table and new QC plots to
    ``restack_dir``. Frames are only excluded here, never moved.

    Typical loop: look at ``diagnostics/frame_quality/*.pdf``, tighten or
    loosen ``best_fraction`` / ``fwhm_max_px`` / ``roundness_max``, run this
    script, compare the stacks.

    Needs ost_photometry >= 0.6.
"""

############################################################################
#                                  Input                                   #
############################################################################
#   Output directory of ``1_reduce_and_stack.py``
output_dir: str = "output/"

#   Where the new stacks go (created if missing)
restack_dir: str = "output/restack/"

#   Target name written to the stack header (``None`` keeps the original)
target_name: str | None = None

############################################################################
#                       Frame quality and selection                        #
############################################################################
#   Same meaning as in ``1_reduce_and_stack.py``; ``None`` = off.
fwhm_max_px: float | None = None
fwhm_max_arcsec: float | None = None
best_fraction: float | None = 0.6
rank_by: str = "fwhm_px"
fwhm_sigma_clip: float | None = None
roundness_max: float | None = 0.3
n_stars_min: int | None = None
background_max: float | None = None
min_frames: int = 5

############################################################################
#                                Stacking                                  #
############################################################################
stack_method: str = "average"
stack_weighting: str = "fwhm"
dtype: str | None = None

############################################################################
#                               Libraries                                  #
############################################################################

import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

from astropy import log

log.setLevel("ERROR")

from ost_photometry import style
from ost_photometry.reduce.frame_selection import (
    FrameSelection,
    group_indices_by_filter,
    mark_selection,
    read_quality_table,
    set_string_column,
    stack_weights,
    summarize_selection,
    write_quality_table,
)
from ost_photometry.reduce.quality import plot_frame_quality
from ost_photometry.reduce.workflow.stack import (
    stack_filter_images,
    stack_meta_for_filter,
)


############################################################################
#                         Helper: build the selection                      #
############################################################################


def build_frame_selection() -> FrameSelection:
    """Translate the parameter block into a ``FrameSelection``."""
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
    return FrameSelection.from_mapping(selection)


############################################################################
#                                  Main                                    #
############################################################################

if __name__ == "__main__":
    start_time = time.time()

    out = Path(output_dir)
    aligned_dir = out / "aligned_lights"
    table_path = out / "frame_quality.ecsv"
    target = Path(restack_dir)

    if not table_path.is_file():
        raise SystemExit(
            f"{style.Bcolors.FAIL}No {table_path}. Run 1_reduce_and_stack.py first."
            f"{style.Bcolors.ENDC}"
        )
    if not aligned_dir.is_dir():
        raise SystemExit(
            f"{style.Bcolors.FAIL}No registered frames in {aligned_dir}. Re-run "
            "1_reduce_and_stack.py with keep_aligned_lights = True."
            f"{style.Bcolors.ENDC}"
        )

    ###
    #   Candidates: frames that exist in aligned_lights
    #
    table = read_quality_table(table_path)
    present = {p.name for p in aligned_dir.iterdir() if p.is_file()}
    keep_rows = np.array([str(f) in present for f in table["file"]], dtype=bool)
    table = table[keep_rows]
    if len(table) == 0:
        raise SystemExit(
            f"{style.Bcolors.FAIL}None of the frames in {table_path} exists in "
            f"{aligned_dir}.{style.Bcolors.ENDC}"
        )
    print(f"   {len(table)} registered frames found in {aligned_dir}")

    ###
    #   New selection and weights (previous rejections are reset)
    #
    selection = build_frame_selection()
    table["rejected"] = np.zeros(len(table), dtype=bool)
    set_string_column(table, "reject_reason", [""] * len(table))
    table["aligned"] = np.ones(len(table), dtype=bool)
    mark_selection(table, selection)
    table["stack_weight"] = stack_weights(table, stack_weighting)
    print(f"   Selection: {selection.describe()}; weighting: {stack_weighting}")
    summarize_selection(table)

    target.mkdir(parents=True, exist_ok=True)
    write_quality_table(table, target / "frame_quality_restack.ecsv")

    ###
    #   Stack per filter
    #
    rejected = np.asarray(table["rejected"], dtype=bool)
    weights_all = np.asarray(table["stack_weight"], dtype=float)
    files_all = [str(f) for f in table["file"]]
    written: list[str] = []
    for filter_, idx in group_indices_by_filter(table).items():
        kept = idx[~rejected[idx]]
        if kept.size == 0:
            print(f"   Filter {filter_}: no frames left after selection, skipped.")
            continue
        paths = [str(aligned_dir / files_all[i]) for i in kept]
        weights = None if stack_weighting == "none" else weights_all[kept]
        meta = stack_meta_for_filter(
            table, filter_, weighting=stack_weighting, stacked_files=paths
        )
        name = stack_filter_images(
            paths,
            stack_method,
            dtype,
            filter_,
            target,
            target_name,
            weights=weights,
            stack_meta=meta,
        )
        written.append(name)
        print(f"   Filter {filter_}: stacked {kept.size} of {idx.size} frames -> {name}")

    plot_frame_quality(table, target, selection=selection, blocking=True)

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print(f"   Stacks:        {', '.join(written) if written else 'none written'}")
    print(f"   Frame quality: {target / 'frame_quality_restack.ecsv'}")
    print(f"   QC plots:      {target / 'diagnostics' / 'frame_quality'}")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
