#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Archive pipeline, optional step 5: colour images of every target from
    the stacks of step 3 (``stacks/``) or step 4 (``restack/``).

    For each target the per-filter stacks (all filters on one grid) are
    mixed into R, G and B, the sky is neutralised, optionally gradients are
    removed and the colours balanced on the stars, then stretched and
    written as 16-bit TIFF plus PNG / JPEG previews.

    Palettes: ``rgb`` (R, G or V, B), ``lrgb`` (plus Clear / L as
    luminance), ``sho`` (SII, Halpha, OIII), ``hoo`` (Halpha, OIII, OIII),
    ``custom`` (free mix), ``auto`` (from the filters of each target).

    Output: ``<output_dir>/color/<target>/<target>_<palette>[_<camera>].tif /
    .png / .jpg`` (camera only when a target has images of several cameras);
    the TIFF description holds the processing log (JSON).

    Needs ost_photometry >= 0.6.
"""

############################################################################
#                                  Input                                   #
############################################################################
#   Output directory of 3_reduce_and_stack.py.
output_dir: str = "output/"

#   "stacks" (step 3, ``<output_dir>/stacks/summary.ecsv``) or "restack"
#   (step 4, ``restack_dir/summary.ecsv``).
source: str = "stacks"
restack_dir: str = "output/restack/"

#   Only these targets (names as in summary.ecsv); ``None`` = all.
targets: list[str] | None = None

#   "auto": the camera combination when it exists, else one image per
#   camera; "combined" or a camera id (e.g. "qhy600m"): only that one.
camera: str = "auto"

############################################################################
#                                 Palette                                  #
############################################################################
#   "auto", "rgb", "lrgb", "sho", "hoo" or "custom".
palette: str = "auto"

#   For palette = "custom": filters and weights per channel, e.g.
custom_mix: dict | None = None
# custom_mix: dict | None = {
#     "R": {"Ha": 0.7, "R": 0.3},
#     "G": {"V": 1.0},
#     "B": {"B": 1.0},
# }
luminance: str | None = None

############################################################################
#                               Processing                                 #
############################################################################
#   Same meaning as in astro_imaging/3_color_composite.py.
neutralize_background: bool = True
crop_to_common_field: bool = True
cosmetic_correction: bool = True
remove_gradient: bool = False
gradient_method: str = "polynomial"
gradient_degree: int = 2
gradient_box: int | None = None
star_color_calibration: bool = False
channel_factors: dict | None = None

############################################################################
#                                 Stretch                                  #
############################################################################
stretch: str = "auto"
linked_stretch: bool | None = None
target_background: float = 0.25
shadows_sigma: float = 2.8
highlight_percentile: float = 99.995
saturation: float = 1.0

############################################################################
#                                  Output                                  #
############################################################################
formats: list[str] = ["tiff", "png", "jpeg"]

############################################################################
#                               Libraries                                  #
############################################################################

import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from astropy import log
from astropy.table import Table

log.setLevel("ERROR")

from ost_photometry import style
from ost_photometry.archive.cache import safe_name
from ost_photometry.imaging import (
    ColorSettings,
    compose,
    stacks_from_summary,
    write_images,
)

############################################################################
#                                  Main                                    #
############################################################################


def color_settings() -> ColorSettings:
    return ColorSettings(
        palette=palette, custom_mix=custom_mix, luminance=luminance,
        neutralize_background=neutralize_background,
        crop_to_common_field=crop_to_common_field,
        cosmetic_correction=cosmetic_correction,
        remove_gradient=remove_gradient, gradient_method=gradient_method,
        gradient_degree=gradient_degree, gradient_box=gradient_box,
        star_color_calibration=star_color_calibration, channel_factors=channel_factors,
        stretch=stretch, linked_stretch=linked_stretch,
        target_background=target_background, shadows_sigma=shadows_sigma,
        highlight_percentile=highlight_percentile, saturation=saturation,
    )


if __name__ == "__main__":
    start_time = time.time()
    out = Path(output_dir).expanduser()
    summary_path = (out / "stacks" if source == "stacks"
                    else Path(restack_dir).expanduser()) / "summary.ecsv"
    if not summary_path.is_file():
        raise SystemExit(f"{style.Bcolors.FAIL}No {summary_path}. Run "
                         f"{'3_reduce_and_stack.py' if source == 'stacks' else '4_restack.py'}"
                         f" first.{style.Bcolors.ENDC}")
    summary = Table.read(summary_path, format="ascii.ecsv")
    wanted = {t.lower().replace(" ", "") for t in targets} if targets else None
    settings = color_settings()

    written = []
    stack_sets = stacks_from_summary(summary, camera=camera)
    per_target: dict[str, int] = {}
    for stack_set in stack_sets:
        per_target[stack_set.target] = per_target.get(stack_set.target, 0) + 1
    for stack_set in stack_sets:
        if wanted and stack_set.target.lower().replace(" ", "") not in wanted:
            continue
        label = f"{stack_set.target} ({stack_set.camera})"
        try:
            result = compose(stack_set.stacks, settings)
        except ValueError as exc:
            print(f"   {label}: skipped - {exc}")
            continue
        # camera in the name only when the target has images of several cameras
        suffix = f"_{stack_set.camera}" if per_target[stack_set.target] > 1 else ""
        base = (out / "color" / safe_name(stack_set.target)
                / f"{safe_name(stack_set.target)}_{result.palette.name}{suffix}")
        files = write_images(result, base, formats=formats)
        written.extend(files)
        factors = ", ".join(f"{k} {v:.2f}" for k, v in result.log["channel_factors"].items())
        print(f"   {label}: {result.palette.name} from {', '.join(sorted(stack_set.stacks))}"
              f"; factors {factors}")
        if "hint" in result.log:
            print(f"      hint: {result.log['hint']}")

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    for path in written:
        print(f"   Image: {path}")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
