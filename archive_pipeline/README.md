# Archive pipeline: fetch, group, reduce, stack

Three scripts that take an object or an observation run from the
[OST data archive](https://polaris.astro.physik.uni-potsdam.de/data_archive/),
work out which calibration frames belong to which science frames, reduce
everything with the right masters and stack every target.

Needs `ost_photometry` ≥ 0.6 and ASTAP (`astap_cli`) for plate solving.

| Script | Role |
|--------|------|
| `1_fetch.py` | Download the science frames of an object or run, all bias / dark / flat candidates of these runs and of neighbouring runs (`calib_window_days`), and metadata of other lights in that window. Science exposures without matching darks get darks from the archive's dark finder (`use_dark_finder`, needs a login). Content-addressed cache, checksums verified. Works also on a local directory tree (`local_directory`). |
| `2_classify_and_group.py` | Frame types from image statistics, electronic setups, targets by sky position, camera orientation (archive WCS or local ASTAP), mount sessions, flat probabilities, reduction units. Writes `calibration_plan.yaml`, `calibration_groups.ecsv` and timeline plots. |
| `3_reduce_and_stack.py` | Masters per calibration group, lights per unit, then per target: quality selection, registration onto one grid, weighted stack per camera and filter over all nights; optional camera combination. |

```bash
export OST_ARCHIVE_USER=...        # or answer the prompt
python 1_fetch.py                  # run_name = "2022-03-08" or object_name = "M57"
python 2_classify_and_group.py
# check grouping_report.txt and diagnostics/calibration_groups/*.pdf,
# correct the 'overrides' block in calibration_plan.yaml if needed, re-run step 2
python 3_reduce_and_stack.py
```

## Concept

Calibration and stacking use different groupings.

| Level | Key | Used for |
|-------|-----|----------|
| Electronic setup | camera, binning, readout mode, gain, offset, temperature (±2 K) | bias and darks; reusable across nights |
| Mount session | camera, telescope, orientation mod 180°, parity, pixel scale | time block without remounting |
| Flat group | camera, binning, filter, mount session | flats |
| Target | cluster of field centres (not the `OBJECT` name) | reference frame, selection, stack |

Calibration runs per **reduction unit** (mount session × electronic setup):
all targets of a night share the same masters. Stacking runs per **target ×
camera × filter** over all sessions and nights in which the target appears.
All filters and cameras of a target are registered onto one grid.

**Why orientation and not flats?** Plate solutions of the test data show the
camera orientation stable to 0.2° within a night (also across meridian
flips) but changing by 8°–66° between nights. Flats reveal a rotation only
above roughly 90°: their sharp dust donuts sit on the sensor window and
rotate with the camera. The dust pattern is still useful as a fingerprint of
"same camera, same dust state".

## Flat probabilities

| Category | Meaning |
|----------|---------|
| certain | flats taken during the session, or p ≥ 0.95 |
| likely | p ≥ 0.7, e.g. dawn flats after the session with matching dust pattern |
| uncertain | p ≥ 0.4, used with a warning |
| rejected | other session or camera in between, different binning, or too far away |

p starts from `exp(-Δt/τ)` (τ = 2 days) and is multiplied by likelihood
factors: dust pattern matches (×3) or differs (×0.2), vignetting differs
(×0.1). Flats before and after a session are combined when both are at
least "likely". Without an applicable flat, `no_flat_policy` decides
(`best_available`, `skip_flat`, `exclude_lights`).

## Overrides (`calibration_plan.yaml`)

```yaml
overrides:
  exclude_frames: ["123456"]                      # archive pk / frame id
  frame_types: {"123457": flat}
  session_breaks: ["123500"]                      # remount the orientation did not show
  merge_sessions: [[S20220308_01_qhy600m, S20220309_01_qhy600m]]
  force_flats: {S20220308_01_qhy600m: {V: [FS_20220308_qhy600m_3x3_V]}}
  merge_targets: [["M31 panel 1", "M31 panel 2"]]
  rename_targets: {T03: "NGC 7000"}
  no_stack_targets: ["field_283.396+33.029"]
```

Edit, then run `2_classify_and_group.py` again (plate solutions are cached).

## Output

| Path | Content |
|------|---------|
| `<workspace>/manifest.ecsv` | all frames with archive and header metadata |
| `<workspace>/calibration_groups.ecsv` | frame type, setup, session, target, flat set, masters per frame |
| `<workspace>/calibration_plan.yaml` | units, masters, targets, flat candidates with probabilities |
| `<workspace>/diagnostics/calibration_groups/timeline_*.pdf` | orientation and frames over the night |
| `<output>/masters/`, `<output>/reduced/<unit>/` | masters and reduced lights (e-/s) |
| `<output>/stacks/<target>/<camera>/combined_filter_<F>.fit` | stacks |
| `<output>/stacks/summary.ecsv` | targets, cameras, filters, frames, exposure |

## Limits

- Spectroscopy frames are excluded (archive classification or spectral
  structure); this pipeline is for imaging.
- Moving solar-system objects: set `target_grouping = "archive_object"`;
  object-tracked stacking is not included.
- Nights without solvable lights give unverified sessions; flat assignment
  then relies on time, dust and camera changes only.
- Mount sessions are not written back to the archive.
