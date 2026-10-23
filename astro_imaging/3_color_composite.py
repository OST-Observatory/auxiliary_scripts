#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Astro imaging, step 3: colour image from the per-filter stacks.

    Reads the stacks ``combined_filter_<F>.fit`` written by
    ``1_reduce_and_stack.py`` (or ``2_restack.py``) - all filters on one
    pixel grid - mixes them into R, G and B, neutralises the sky, optionally
    removes gradients and balances the colours on the stars, stretches and
    writes a 16-bit TIFF (for further processing in Siril, PixInsight,
    GIMP) plus PNG / JPEG previews.

    Palettes: ``rgb`` (R, G or V, B), ``lrgb`` (plus Clear / L as
    luminance), ``sho`` (SII, Halpha, OIII), ``hoo`` (Halpha, OIII, OIII),
    ``custom`` (free mix, e.g. Halpha into red), ``auto`` (from the filters).

    Output in ``output_dir``: ``<name>_<palette>.tif / .png / .jpg``; the
    TIFF description holds the processing log (JSON).

    Needs ost_photometry >= 0.6.
"""

############################################################################
#                                  Input                                   #
############################################################################
#   Directory with the stacks: ``output/`` of 1_reduce_and_stack.py or
#   ``output/restack/`` of 2_restack.py.
stacks_dir: str = "output/"

#   Where the images go (``None`` = ``<stacks_dir>/color/``).
output_dir: str | None = None

#   File name stem (``None`` = OBJECT keyword of the stacks).
name: str | None = None

############################################################################
#                                 Palette                                  #
############################################################################
#   "auto", "rgb", "lrgb", "sho", "hoo" or "custom".
palette: str = "auto"

#   For palette = "custom": filters and weights per channel. Filter names as
#   in the headers or short (Ha, OIII, SII, L for Clear), e.g. Halpha into
#   the red channel of an RGB image:
custom_mix: dict | None = None
# custom_mix: dict | None = {
#     "R": {"Ha": 0.7, "R": 0.3},
#     "G": {"V": 1.0},
#     "B": {"B": 1.0},
# }
#   Luminance filter for palette = "custom" (e.g. "Clear"); ``None`` = none.
luminance: str | None = None

############################################################################
#                               Processing                                 #
############################################################################
#   Set the sky of every channel to the same (grey) level.
neutralize_background: bool = True

#   Cut the image to the field all filters cover (no coloured edges).
crop_to_common_field: bool = True

#   Replace cold / dead pixels (black dots after stretching).
cosmetic_correction: bool = True

#   Remove large-scale gradients (moon, light pollution, flat residuals):
#   "polynomial" (smooth 2-D polynomial of ``gradient_degree``, robust next
#   to large galaxies / nebulae) or "background2d" (follows smaller
#   structures; ``gradient_box`` in pixels, default 1/8 of the shorter side).
#   Rings and blotches from flats that do not quite match: "background2d"
#   with ``gradient_box`` about 1/12 of the shorter side (e.g. 177 for
#   3192 x 2129). Stars and objects are masked; under a large object the
#   sky is interpolated, so a faint halo can remain there. Better flats
#   (same night, same camera position) remove the cause.
remove_gradient: bool = False
gradient_method: str = "polynomial"
gradient_degree: int = 2
gradient_box: int | None = None

#   White balance so that the average star is neutral (broadband palettes).
star_color_calibration: bool = False

#   Or fixed channel factors, e.g. {"R": 0.8, "G": 1.0, "B": 1.4}
#   (``None`` = 1, or the star calibration).
channel_factors: dict | None = None

############################################################################
#                                 Stretch                                  #
############################################################################
#   "auto" (midtone transfer: sky at ``target_background``), "asinh",
#   "linear".
stretch: str = "auto"

#   One stretch for all channels (keeps broadband colours) or one per
#   channel (narrowband). ``None`` = by palette.
linked_stretch: bool | None = None

#   Brightness of the sky after the stretch (0.1 dark ... 0.3 bright).
target_background: float = 0.25

#   Black point at sky - shadows_sigma x noise; white point percentile.
shadows_sigma: float = 2.8
highlight_percentile: float = 99.995

#   Colour saturation (1 = unchanged).
saturation: float = 1.0

############################################################################
#                                  Output                                  #
############################################################################
#   Any of "tiff" (16 bit), "png", "jpeg".
formats: list[str] = ["tiff", "png", "jpeg"]

############################################################################
#                               Libraries                                  #
############################################################################

import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from astropy import log
from astropy.io import fits

log.setLevel("ERROR")

from ost_photometry import style
from ost_photometry.imaging import (
    ColorSettings,
    compose,
    stacks_from_directory,
    write_images,
)
from ost_photometry.archive.cache import safe_name

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
    source = Path(stacks_dir).expanduser()
    stacks = stacks_from_directory(source)
    if not stacks:
        raise SystemExit(f"{style.Bcolors.FAIL}No combined_filter_*.fit in {source}. "
                         f"Run 1_reduce_and_stack.py first.{style.Bcolors.ENDC}")
    print(f"   Stacks: {', '.join(sorted(stacks))}")
    result = compose(stacks, color_settings())
    stem = name or str(fits.getheader(next(iter(stacks.values()))).get("OBJECT", "") or "image")
    target_dir = Path(output_dir).expanduser() if output_dir else source / "color"
    files = write_images(result, target_dir / f"{safe_name(stem)}_{result.palette.name}",
                         formats=formats)

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    mix = "; ".join(f"{ch}=" + "+".join(f"{w:g}x{f}" for f, w in weights.items())
                    for ch, weights in result.palette.mix.items())
    print(f"   Palette:  {result.palette.name} ({mix}"
          + (f"; L={result.palette.luminance}" if result.palette.luminance else "") + ")")
    print(f"   Factors:  {result.log['channel_factors']}"
          + (f" from {result.log['calibration_stars']} stars"
             if "calibration_stars" in result.log else ""))
    if "hint" in result.log:
        print(f"   Hint:     {result.log['hint']}")
    for path in files:
        print(f"   Image:    {path}")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
