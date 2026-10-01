#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Archive pipeline, step 2 of 3: classify and group the fetched frames.

      * frame types from image statistics (header IMAGETYP is often wrong),
        spectroscopy frames are recognised and excluded,
      * electronic setups (camera, binning, readout, gain, offset,
        temperature) for bias / darks,
      * targets by sky position (several targets per night are normal),
      * camera orientation of sampled lights (archive plate solution when
        available, otherwise local ASTAP) and mount sessions,
      * flat sets with the probability that they belong to a session,
      * reduction units (session x electronic setup) and their masters.

    Output in ``workspace``:
      * ``calibration_groups.ecsv``  one row per frame with all group ids
      * ``calibration_plan.yaml``    units, masters, targets, probabilities;
                                     edit the ``overrides`` block and run
                                     this script again to apply corrections
      * ``diagnostics/calibration_groups/timeline_<night>_<telescope>.pdf``
      * ``orientation_cache.ecsv``   plate solutions (re-runs are fast)

    Needs ost_photometry >= 0.6 and ASTAP (astap_cli) for plate solving.
"""

############################################################################
#                                  Input                                   #
############################################################################
#   Working directory of ``1_fetch.py``.
workspace: str = "archive_workspace/"

############################################################################
#                Additional options: only edit if necessary                #
############################################################################
#   Orientation change (deg, modulo 180) that starts a new mount session.
#   Within one night the orientation is stable to about 0.2 deg.
pa_tolerance: float = 0.5

#   Plate-solve one light frame per interval (minutes) plus block ends,
#   gaps and target changes; changes are bisected automatically.
solve_interval_minutes: float = 30.0

#   Seconds before an ASTAP run is abandoned.
solve_timeout: float = 120.0

#   Two frames belong to the same target if their field centres are closer
#   than this fraction of the field height. ``target_grouping`` = "position"
#   (default), "name" or "archive_object" (moving solar-system objects).
target_overlap_fraction: float = 0.5
target_grouping: str = "position"

#   Flat probability: time constant of the decay (days). Flats taken
#   ``tau_days`` away from a session start with p = exp(-1) before the dust
#   and vignetting checks.
tau_days: float = 2.0

#   Probability limits of the categories certain / likely / uncertain.
p_certain: float = 0.95
p_likely: float = 0.7
p_uncertain: float = 0.4

#   What to do without an applicable flat: "best_available" (use the best
#   candidate and warn), "skip_flat" (reduce without flat field),
#   "exclude_lights" (do not reduce those lights).
no_flat_policy: str = "best_available"

#   Bias / darks: temperature tolerance (K), exposure tolerance (s), and
#   the search window (days).
temp_tolerance: float = 2.0
dark_exptime_tolerance: float = 0.5
calibration_window_days: float = 30.0

#   Worker processes for the image statistics (``None`` = half the CPUs).
n_cores_multiprocessing: int | None = None

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
from ost_photometry.archive import read_manifest
from ost_photometry.reduce.grouping.plan import (
    PlanSettings,
    build_calibration_plan,
    read_overrides,
    write_frames,
    write_plan,
)
from ost_photometry.reduce.grouping.plots import plot_night_timelines

############################################################################
#                                  Main                                    #
############################################################################


def _progress(i: int, n: int, name: str) -> None:
    print(f"   [{i}/{n}] plate solving {name}", flush=True)


if __name__ == "__main__":
    start_time = time.time()
    work = Path(workspace).expanduser()
    manifest_path = work / "manifest.ecsv"
    if not manifest_path.is_file():
        raise SystemExit(f"{style.Bcolors.FAIL}No {manifest_path}; run 1_fetch.py first."
                         f"{style.Bcolors.ENDC}")
    plan_path = work / "calibration_plan.yaml"
    overrides = read_overrides(plan_path)
    settings = PlanSettings(
        temp_tolerance=temp_tolerance,
        pa_tolerance=pa_tolerance,
        solve_interval_minutes=solve_interval_minutes,
        solve_timeout=solve_timeout,
        target_overlap_fraction=target_overlap_fraction,
        target_grouping=target_grouping,
        dark_exptime_tolerance=dark_exptime_tolerance,
        calibration_window_days=calibration_window_days,
        no_flat_policy=no_flat_policy,
        n_cores_multiprocessing=n_cores_multiprocessing,
        flat={"tau_days": tau_days, "p_certain": p_certain, "p_likely": p_likely,
              "p_uncertain": p_uncertain},
    )
    plan = build_calibration_plan(
        read_manifest(manifest_path),
        work_dir=work,
        settings=settings,
        overrides=overrides,
        progress=_progress,
        log=lambda text: print("   " + text, flush=True),
    )
    write_plan(plan, plan_path)
    write_frames(plan, work / "calibration_groups.ecsv")
    plots = plot_night_timelines(plan, work)
    (work / "grouping_report.txt").write_text("\n".join(plan.report) + "\n")
    for line in plan.report:
        print("   " + line)
    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print(f"   Plan:      {plan_path}  (edit 'overrides', then run this script again)")
    print(f"   Frames:    {work / 'calibration_groups.ecsv'}")
    print(f"   Timelines: {len(plots)} PDF(s) in {work / 'diagnostics' / 'calibration_groups'}")
    print("   Next: 3_reduce_and_stack.py")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
