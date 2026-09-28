#!/usr/bin/env python3
"""Measure frame quality (FWHM, roundness, stars, sky) for a directory of FITS.

Standalone front end to ``ost_photometry.reduce.quality``: works on reduced
lights, on ``aligned_lights/`` of a reduction, or on a sequence registered by
Siril. Without selection flags it only measures and writes the table; with
flags it marks rejections (and can move the files), computes stack weights,
writes header keywords, and draws the per-filter QC plot.

Usage::

    python frame_quality.py output/aligned_lights
    python frame_quality.py lights/ --best-fraction 0.8 --roundness-max 0.3 --plot
    python frame_quality.py lights/ --fwhm-max 3.5 --move-rejected lights/rejected
    python frame_quality.py siril_seq/ --imagetyp LIGHT --weighting fwhm --write-headers

Needs ost_photometry >= 0.6 (installed, or checked out next to this repo).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Defaults (also editable here; every one is an argparse default)
# ---------------------------------------------------------------------------
PATH_INPUT = "?"
PATH_OUTPUT_TABLE: str | None = None  # default: <input>/frame_quality.ecsv
IMAGETYP: str | None = None  # e.g. "LIGHT"; None = every FITS in the directory
FWHM_MAX: float | None = None
FWHM_MAX_ARCSEC: float | None = None
BEST_FRACTION: float | None = None
RANK_BY = "fwhm_px"
SIGMA_CLIP: float | None = None
ROUNDNESS_MAX: float | None = None
N_STARS_MIN: int | None = None
BACKGROUND_MAX: float | None = None
MASKED_FRACTION_MAX: float | None = None
MIN_FRAMES = 1
WEIGHTING = "none"
N_CORES: int | None = None

# Prefer a source checkout next to this repository, else the installed package.
_PROJECTS = Path(__file__).resolve().parents[2]
_PKG_SRC = _PROJECTS / "ost_photometry_package" / "src"
if _PKG_SRC.is_dir() and str(_PKG_SRC) not in sys.path:
    sys.path.insert(0, str(_PKG_SRC))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("input_dir", nargs="?", default=PATH_INPUT, help="Directory with FITS frames")
    p.add_argument("--out", default=PATH_OUTPUT_TABLE, help="ECSV table path")
    p.add_argument("--imagetyp", default=IMAGETYP, help="Restrict to this IMAGETYP")
    sel = p.add_argument_group("selection (all set criteria must pass, per filter)")
    sel.add_argument("--fwhm-max", type=float, default=FWHM_MAX, help="max FWHM [pixel]")
    sel.add_argument(
        "--fwhm-max-arcsec", type=float, default=FWHM_MAX_ARCSEC, help="max FWHM [arcsec]"
    )
    sel.add_argument(
        "--best-fraction", type=float, default=BEST_FRACTION, help="keep the best fraction (0-1]"
    )
    sel.add_argument(
        "--rank-by", default=RANK_BY, choices=["fwhm_px", "fwhm_weighted"],
        help="ranking for --best-fraction",
    )
    sel.add_argument(
        "--sigma-clip", type=float, default=SIGMA_CLIP, help="reject FWHM > median + k MAD"
    )
    sel.add_argument("--roundness-max", type=float, default=ROUNDNESS_MAX)
    sel.add_argument("--n-stars-min", type=int, default=N_STARS_MIN)
    sel.add_argument("--background-max", type=float, default=BACKGROUND_MAX)
    sel.add_argument("--masked-fraction-max", type=float, default=MASKED_FRACTION_MAX)
    sel.add_argument("--min-frames", type=int, default=MIN_FRAMES, help="floor per filter")
    p.add_argument(
        "--weighting", default=WEIGHTING, choices=["none", "fwhm", "n_stars", "noise"],
        help="stack weights written to the table (and FRMWGHT with --write-headers)",
    )
    p.add_argument(
        "--write-headers", action="store_true",
        help="write FWHM / NSTARS / ... and FRMWGHT into the frame headers",
    )
    p.add_argument(
        "--move-rejected", metavar="DIR", default=None,
        help="move rejected frames to DIR (header gets QCREJ / QCREASON)",
    )
    p.add_argument(
        "--plot", nargs="?", const="__input__", default=None, metavar="DIR",
        help="write diagnostics/frame_quality/*.pdf below DIR (default: input directory)",
    )
    p.add_argument("--ncores", type=int, default=N_CORES, help="worker processes (1 = serial)")
    return p.parse_args(argv)


def build_selection(args: argparse.Namespace):
    from ost_photometry.reduce.frame_selection import FrameSelection

    if args.fwhm_max is not None and args.fwhm_max_arcsec is not None:
        raise SystemExit("Use either --fwhm-max or --fwhm-max-arcsec, not both.")
    mapping: dict = {"min_frames": args.min_frames, "rank_by": args.rank_by}
    if args.fwhm_max is not None:
        mapping.update(fwhm_max=args.fwhm_max, fwhm_unit="px")
    elif args.fwhm_max_arcsec is not None:
        mapping.update(fwhm_max=args.fwhm_max_arcsec, fwhm_unit="arcsec")
    for key, value in (
        ("best_fraction", args.best_fraction),
        ("fwhm_sigma_clip", args.sigma_clip),
        ("roundness_max", args.roundness_max),
        ("n_stars_min", args.n_stars_min),
        ("background_max", args.background_max),
        ("masked_fraction_max", args.masked_fraction_max),
    ):
        if value is not None:
            mapping[key] = value
    return FrameSelection.from_mapping(mapping)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.input_dir in (None, "?"):
        print("Give the input directory (positional) or set PATH_INPUT in the script.")
        return 2
    input_dir = Path(args.input_dir).expanduser()
    if not input_dir.is_dir():
        print(f"Not a directory: {input_dir}")
        return 2

    from ost_photometry.reduce.frame_selection import (
        mark_selection,
        stack_weights,
        summarize_selection,
        write_quality_table,
    )
    from ost_photometry.reduce.quality import (
        measure_directory_quality,
        move_rejected_frames,
        plot_frame_quality,
        write_quality_to_headers,
        write_weights_to_headers,
    )

    selection = build_selection(args)
    print(f"Measuring frame quality in {input_dir} ...", flush=True)
    table = measure_directory_quality(
        input_dir,
        image_type_list=[args.imagetyp] if args.imagetyp else None,
        n_cores_multiprocessing=args.ncores,
    )
    if len(table) == 0:
        print("No FITS frames found (check --imagetyp).")
        return 2

    if selection.is_active():
        print(f"Selection: {selection.describe()}")
    mark_selection(table, selection)
    table["stack_weight"] = stack_weights(table, args.weighting)
    summarize_selection(table)

    table_path = Path(args.out) if args.out else input_dir / "frame_quality.ecsv"
    write_quality_table(table, table_path)
    print(f"Table: {table_path}")

    if args.write_headers:
        n_written = write_quality_to_headers(table, input_dir)
        n_weights = write_weights_to_headers(table, input_dir)
        print(f"Headers updated: {n_written} frames ({n_weights} with FRMWGHT)")

    if args.plot is not None:
        plot_root = input_dir if args.plot == "__input__" else Path(args.plot).expanduser()
        paths = plot_frame_quality(table, plot_root, selection=selection, blocking=True)
        for path in paths:
            print(f"Plot: {path}")

    if args.move_rejected:
        moved = move_rejected_frames(table, input_dir, Path(args.move_rejected).expanduser())
        print(f"Moved {len(moved)} rejected frame(s) to {args.move_rejected}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
