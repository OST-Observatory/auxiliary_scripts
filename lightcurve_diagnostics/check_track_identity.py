#!/usr/bin/env python3
"""Identity vs calibration diagnostics for a finished C7 pipeline run.

Reads ``<output_dir>/tables/calibrated_magnitudes*.ecsv`` (or the instrumental
``extracted_magnitudes*.ecsv`` fallback) and answers two questions per star:

1. Is the track one star?  On a registered series (``aa_true`` / ``wcs`` with
   ``shift_all=True``) the same ``id`` must sit at the same pixel in every
   epoch.  A large x/y spread means several stars were mixed into one id.

2. Is the jump instrumental or introduced by the calibration?  Compares the
   epoch-to-epoch scatter of the instrumental magnitude ``mag_<F>`` with that
   of ``mag_cal_<F>``.  A star that is stable in ``mag_<F>`` but jumps in
   ``mag_cal_<F>`` points at the per-epoch zero point / colour term.

Per epoch it also prints the number of calibrators used and the ZP proxy
``median(mag_std - mag_inst)`` over those calibrators.

Usage::

    python check_track_identity.py output/ --filter V --ids 392 --top 3
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
from astropy.table import Table


def _load_table(output_dir: str) -> Table:
    tdir = os.path.join(output_dir, "tables")
    for stem in ("calibrated_magnitudes", "extracted_magnitudes"):
        hits = sorted(glob.glob(os.path.join(tdir, f"{stem}*.ecsv")))
        if hits:
            print(f"Reading {hits[0]}")
            return Table.read(hits[0], format="ascii.ecsv")
    sys.exit(f"No calibrated_magnitudes*.ecsv / extracted_magnitudes*.ecsv in {tdir}")


def _col(tbl: Table, *names: str) -> np.ndarray | None:
    for n in names:
        if n in tbl.colnames:
            return np.asarray(tbl[n], dtype=float)
    return None


def _epoch_jd(tbl: Table) -> np.ndarray | None:
    return _col(tbl, "observation_jd", "jd")


def _nanmedian(a: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else np.nan


def _nanstd(a: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    return float(np.std(a)) if a.size > 1 else np.nan


def _fmt(v: float, w: int = 7, p: int = 3) -> str:
    return f"{v:{w}.{p}f}" if np.isfinite(v) else " " * (w - 3) + "nan"


def report_epochs(tbl: Table, filt: str) -> None:
    eids = np.asarray(tbl["epoch_id"]).astype(str)
    inst = _col(tbl, f"mag_{filt}", f"mag_inst_{filt}")
    std = _col(tbl, f"mag_std_{filt}")
    cal = _col(tbl, f"mag_cal_{filt}")
    used = (
        np.asarray(tbl[f"is_calibrator_{filt}"], dtype=bool)
        if f"is_calibrator_{filt}" in tbl.colnames
        else None
    )
    jd = _epoch_jd(tbl)
    print(f"\n=== Epoch summary, filter {filt} ===")
    print("epoch        jd            n_rows  n_cal  ZP_proxy  ZP_std   cal-inst_med  cal-inst_std")
    for eid in sorted(set(eids)):
        m = eids == eid
        n_rows = int(np.sum(m))
        jd_s = f"{np.nanmedian(jd[m]):.5f}" if jd is not None else "?"
        if inst is None or std is None:
            print(f"{eid:12s} {jd_s:>13s} {n_rows:6d}   (no instrumental / catalog columns)")
            continue
        sel = m & np.isfinite(inst) & np.isfinite(std)
        if used is not None:
            sel &= used
        n_cal = int(np.sum(sel))
        zp = std[sel] - inst[sel]
        zp_med = float(np.median(zp)) if n_cal else np.nan
        zp_std = float(np.std(zp)) if n_cal > 1 else np.nan
        if cal is not None:
            d = cal[sel] - inst[sel]
            d_med = float(np.median(d)) if n_cal else np.nan
            d_std = float(np.std(d)) if n_cal > 1 else np.nan
        else:
            d_med = d_std = np.nan
        print(
            f"{eid:12s} {jd_s:>13s} {n_rows:6d} {n_cal:6d}  {_fmt(zp_med)}  {_fmt(zp_std)}"
            f"     {_fmt(d_med)}      {_fmt(d_std)}"
        )
    print(
        "\nZP_std is the scatter of (catalog - instrumental) over the used calibrators;"
        "\ncal-inst_std > 0.05 means the colour term moves stars by different amounts."
    )


def report_ids(tbl: Table, filt: str, ids: list[int], verbose: bool) -> None:
    all_ids = np.asarray(tbl["id"], dtype=int)
    eids = np.asarray(tbl["epoch_id"]).astype(str)
    x = _col(tbl, "x")
    y = _col(tbl, "y")
    inst = _col(tbl, f"mag_{filt}", f"mag_inst_{filt}")
    cal = _col(tbl, f"mag_cal_{filt}")
    std = _col(tbl, f"mag_std_{filt}")
    jd = _epoch_jd(tbl)
    am = _col(tbl, f"airmass_{filt}", "airmass")
    print(f"\n=== Per-star check, filter {filt} ===")
    print("id      n_ep  dx_rms  dy_rms  max_off   inst_rms  cal_rms   cat_mag  cal_med")
    for sid in ids:
        m = all_ids == int(sid)
        if not np.any(m):
            print(f"{sid:<7d} not in table")
            continue
        n_ep = int(np.sum(m))
        if x is not None and y is not None:
            dx = x[m] - np.median(x[m])
            dy = y[m] - np.median(y[m])
            dx_rms = float(np.sqrt(np.mean(dx**2)))
            dy_rms = float(np.sqrt(np.mean(dy**2)))
            max_off = float(np.max(np.hypot(dx, dy)))
        else:
            dx_rms = dy_rms = max_off = np.nan
        inst_rms = _nanstd(inst[m]) if inst is not None else np.nan
        cal_rms = _nanstd(cal[m]) if cal is not None else np.nan
        cat = _nanmedian(std[m]) if std is not None else np.nan
        cal_med = _nanmedian(cal[m]) if cal is not None else np.nan
        flag = ""
        if np.isfinite(max_off) and max_off > 3.0:
            flag += "  <-- position jumps: several stars in one id"
        if np.isfinite(inst_rms) and np.isfinite(cal_rms) and cal_rms > 2.0 * inst_rms + 0.02:
            flag += "  <-- calibration adds scatter (ZP / colour term)"
        print(
            f"{sid:<7d} {n_ep:4d}  {_fmt(dx_rms, 6, 2)}  {_fmt(dy_rms, 6, 2)}  {_fmt(max_off, 7, 2)}"
            f"   {_fmt(inst_rms)}   {_fmt(cal_rms)}   {_fmt(cat)}  {_fmt(cal_med)}{flag}"
        )
        if verbose:
            order = np.argsort(jd[m]) if jd is not None else np.arange(n_ep)
            print("    epoch        jd          x        y      airmass  inst     cal      cat")
            for k in order:
                row = np.flatnonzero(m)[k]
                print(
                    f"    {eids[row]:12s} {_fmt(jd[row] if jd is not None else np.nan, 11, 4)}"
                    f" {_fmt(x[row] if x is not None else np.nan, 8, 1)}"
                    f" {_fmt(y[row] if y is not None else np.nan, 8, 1)}"
                    f"  {_fmt(am[row] if am is not None else np.nan, 6, 2)}"
                    f"  {_fmt(inst[row] if inst is not None else np.nan)}"
                    f"  {_fmt(cal[row] if cal is not None else np.nan)}"
                    f"  {_fmt(std[row] if std is not None else np.nan)}"
                )


def ooi_and_top_ids(output_dir: str, filt: str, top: int) -> tuple[list[int], list[int]]:
    tdir = os.path.join(output_dir, "tables")
    ooi: list[int] = []
    lc_path = os.path.join(tdir, "light_curves.ecsv")
    if os.path.exists(lc_path):
        lc = Table.read(lc_path, format="ascii.ecsv")
        names = np.asarray(lc["object_name"]).astype(str)
        ids = np.asarray(lc["id"], dtype=int)
        ooi = sorted({int(i) for i, n in zip(ids, names, strict=True) if n.strip()})
    worst: list[int] = []
    st_path = os.path.join(tdir, "calibrator_variability_stats.ecsv")
    if os.path.exists(st_path) and top > 0:
        st = Table.read(st_path, format="ascii.ecsv")
        st = st[np.asarray(st["filter"]).astype(str) == filt]
        if len(st):
            order = np.argsort(-np.asarray(st["excess_rms"], dtype=float))
            worst = [int(v) for v in np.asarray(st["id"])[order][:top]]
    return ooi, worst


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output_dir")
    ap.add_argument("--filter", default="V")
    ap.add_argument("--ids", type=int, nargs="*", default=[], help="extra star ids to inspect")
    ap.add_argument("--top", type=int, default=3, help="worst calibrators by excess RMS")
    ap.add_argument("--all", action="store_true", help="scan every id for position jumps")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-epoch rows for each id")
    args = ap.parse_args()

    tbl = _load_table(args.output_dir)
    report_epochs(tbl, args.filter)

    ooi, worst = ooi_and_top_ids(args.output_dir, args.filter, args.top)
    ids = list(dict.fromkeys(ooi + worst + args.ids))
    if ooi:
        print(f"\nObjects of interest (from light_curves.ecsv): {ooi}")
    if worst:
        print(f"Worst calibrators by excess RMS: {worst}")
    if ids:
        report_ids(tbl, args.filter, ids, args.verbose)

    if args.all:
        x = _col(tbl, "x")
        y = _col(tbl, "y")
        if x is None or y is None:
            print("\nNo x/y columns; cannot scan positions.")
            return
        all_ids = np.asarray(tbl["id"], dtype=int)
        bad: list[tuple[int, float]] = []
        for sid in np.unique(all_ids):
            m = all_ids == sid
            if np.sum(m) < 3:
                continue
            off = np.hypot(x[m] - np.median(x[m]), y[m] - np.median(y[m]))
            if np.max(off) > 3.0:
                bad.append((int(sid), float(np.max(off))))
        bad.sort(key=lambda t: -t[1])
        print(
            f"\nPosition scan: {len(bad)} of {np.unique(all_ids).size} ids move by > 3 px "
            "between epochs"
        )
        for sid, off in bad[:20]:
            print(f"   id {sid}: max offset {off:.1f} px")
        if len(bad) > 0.05 * np.unique(all_ids).size:
            print(
                "   Many ids move: the frames are NOT on one pixel grid, or intra-filter"
                " matching ran on the sky with noisy per-frame WCS. Check the log line"
                " 'Intra-filter matching in pixel coordinates'."
            )


if __name__ == "__main__":
    main()
