# Astro imaging: quality-selected, weighted stacking

Scripts that turn raw frames into the best possible **linear** stack per filter,
all filters on **one pixel grid**, as the starting point for colour
compositing and stretching in Siril, PixInsight, GIMP, … They replace the
Siril sequence steps "register → filter by FWHM / roundness / best X % →
weighted stack" with `ost_photometry.reduce`. Nothing here stretches or
combines colours.

Needs `ost_photometry` ≥ 0.6 (frame quality, weighted stacking).

| Script | Role |
|--------|------|
| `1_reduce_and_stack.py` | Bias/dark/flat, cosmic rays, per-frame quality, selection, registration of every filter onto the sharpest frame (`shift_all=True`), weighted sigma-clipped stack per filter, WCS on the stacks, QC plots. Parameter block at the top, run with `python 1_reduce_and_stack.py`. |
| `2_restack.py` | Re-run only selection + stacking on the registered frames kept by step 1 (`keep_aligned_lights = True`). The Siril loop "change the threshold, restack" without re-reducing. |
| `frame_quality.py` | CLI: measure any directory of FITS (reduced lights, `aligned_lights/`, a Siril-registered sequence), optionally select, weight, write headers, move rejected frames, plot. |

## Usage

```bash
# 1. edit the parameter block (paths, best_fraction, roundness_max, weighting), then
python 1_reduce_and_stack.py

# 2. look at output/diagnostics/frame_quality/frame_quality_<F>.pdf, adjust the
#    selection in 2_restack.py and stack again from output/aligned_lights/
python 2_restack.py

# standalone quality table for any FITS directory
python frame_quality.py output/aligned_lights --best-fraction 0.7 --roundness-max 0.3 --plot
python frame_quality.py siril_registered/ --imagetyp LIGHT --fwhm-max 3.5 --move-rejected siril_registered/rejected
```

## Siril → parameters

| Siril | Here |
|-------|------|
| Filter by FWHM ≤ X | `fwhm_max_px = X` (or `fwhm_max_arcsec`) |
| Filter by wFWHM | `rank_by = "fwhm_weighted"` with `best_fraction` |
| Filter by roundness | `roundness_max = 0.3` (IRAF roundness, 0 = round) |
| Keep the best X % | `best_fraction = X/100` |
| k-sigma filtering | `fwhm_sigma_clip = k` |
| Filter by number of stars / background | `n_stars_min`, `background_max` |
| Weighted stacking (FWHM / number of stars / noise) | `stack_weighting = "fwhm" | "n_stars" | "noise"` |
| Reference image = best frame | `reference_image_selection = "best_fwhm"` (default); `reference_image_index` forces one |
| Rejection sigma clipping in the stack | always on (5σ low / high, MAD) |

Selection criteria are combined with AND and evaluated per filter;
`min_frames` is a floor so a filter never ends up with fewer frames than that.
Frames without any stars (clouds) are rejected as soon as a criterion is set.

## Output (`output/`)

| Path | Content |
|------|---------|
| `combined_filter_<F>.fit` | Linear stack per filter (e-/s, float, sigma-clipped weighted average), same shape and pixel grid for all filters. Header: `N-IMAGES`, `WEIGHTNG`, `NFRAMES0`, `NREJECT`, `NALIGNFL`, `FWHMMED`, `FWHMMAX`, WCS. |
| `frame_quality.ecsv` | One row per light frame: `fwhm_px`, `fwhm_arcsec`, `fwhm_weighted`, `roundness`, `n_stars`, `background`, `background_rms`, `masked_fraction`, `status`, `rejected`, `reject_reason`, `is_reference`, `aligned`, `align_note`, `stack_weight`. |
| `rejected_lights/` | Rejected reduced frames (never deleted), header `QCREJ` / `QCREASON`. |
| `aligned_lights/` | Registered frames (only with `keep_aligned_lights = True`); input for `2_restack.py` or Siril drizzle. Header `FRMWGHT`, `FWHM`, `NSTARS`, … |
| `diagnostics/frame_quality/frame_quality_<F>.pdf` | FWHM, roundness, stars, sky per frame; rejected = red ×, reference = star, not aligned = hollow. |
| `restack/` | Stacks, `frame_quality_restack.ecsv` and plots from `2_restack.py`. |

## Not included (on purpose)

Colour combination (RGB / LRGB / narrowband palettes), background gradient
removal, stretching and TIFF / PNG export. Feed the `combined_filter_*.fit`
files to your post-processing tool of choice.
