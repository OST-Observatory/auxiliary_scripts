# Compare registration methods

Run every shift method implemented in `ost_photometry.reduce.registration` on
**two frames** and plot how they disagree.

Methods (pipeline names):

| Name | What it does |
|------|----------------|
| `skimage` | `phase_cross_correlation` (default in `reduce_main`) |
| `own` | in-house FFT phase correlation |
| `aa` | astroalign, translation only |
| `aa_true` | astroalign similarity (shift + rotation + scale) |
| `flow` | skimage optical flow TV-L1 (dense field) |
| `wcs` | reproject the second frame onto the **reference WCS** (needs a celestial WCS on both; solved with ASTAP in the pipeline if missing) |

`wcs` measures sky alignment, not a blind pixel shift. Frames must carry a real
solution (or `--self-test` attaches a TAN WCS consistent with the synthetic
shift). If the headers share one stale WCS while the pixels have moved, `wcs`
will look aligned on the sky and the residual will stay large — that is the
point of comparing it to `aa_true`.

## Usage

Edit the config block at the top of `compare_registration_methods.py`, or:

```bash
# Two reduced (or unreduced) FITS files
python compare_registration_methods.py \
  --ref /data/night1/light_000.fits \
  --other /data/night1/light_007.fits \
  --out /tmp/reg_compare

# Two directories: pick the n-th FITS in each (sorted)
python compare_registration_methods.py \
  --ref /data/night1/light/ \
  --other /data/night2/light/ \
  --ref-index 0 --other-index 0

# Two raw datasets (bias/dark/flat/light as expected by reduce_main)
python compare_registration_methods.py \
  --ref /data/night1/raw/ \
  --other /data/night2/raw/ \
  --reduce \
  --out /tmp/reg_compare

# Synthetic check (no FITS needed)
python compare_registration_methods.py --self-test --out /tmp/reg_selftest
```

Large chips: `--downsample 4` bins before registration; reported shifts are
scaled back to original pixels. Optical flow is the slowest method.

## Output

| File | Content |
|------|---------|
| `registration_shift_summary.pdf` | Δx/Δy bars, rotation/scale, residual RMS |
| `registration_shift_differences.pdf` | pairwise \|Δshift\| between methods |
| `registration_residuals.pdf` | residual maps side by side |
| `registration_<method>.pdf` | reference / other / aligned / residual |
| `registration_flow_field.pdf` | dense flow (TV-L1 only) |
| `registration_shifts.csv` | numbers for the same |
