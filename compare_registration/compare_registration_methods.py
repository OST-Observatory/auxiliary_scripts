#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""Compare all ``ost_photometry.reduce`` registration methods on two frames.

Point this at two FITS files, two directories of reduced lights, or two raw
datasets (optional calibration via ``reduce_main``). Each implemented shift
method is run and the displacements plus residual images are plotted.
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.visualization import ImageNormalize, ZScaleInterval
from scipy.ndimage import shift as ndi_shift

# ---------------------------------------------------------------------------
# Configuration: edit here, or override on the command line
# ---------------------------------------------------------------------------

# Reference frame (FITS file) or directory of FITS / raw dataset root
PATH_REFERENCE = "?"

# Second frame (shifted relative to the reference)
PATH_OTHER = "?"

# Output directory for PDF/CSV
PATH_OUTPUT = "registration_compare"

# If PATH_* is a directory of FITS, which file to use (0-based, sorted)
FILE_INDEX_REFERENCE = 0
FILE_INDEX_OTHER = 0

# If True, run ``reduce_main`` on each dataset first and pick a light frame
# (each path must be a raw folder as expected by reduce_main).
REDUCE_FIRST = False

# Binning applied before registration (1 = full resolution). Shifts are
# reported in *original* pixels.
DOWNSAMPLE = 1

# Methods to run (subset of the pipeline names)
METHODS = ["skimage", "own", "aa", "aa_true", "flow"]

PLOT_FORMAT = "pdf"

# ---------------------------------------------------------------------------

_PROJECTS = Path(__file__).resolve().parents[2]
_PKG_SRC = _PROJECTS / "ost_photometry_package" / "src"
if _PKG_SRC.is_dir() and str(_PKG_SRC) not in sys.path:
    sys.path.insert(0, str(_PKG_SRC))

FITS_SUFFIXES = {".fit", ".fits", ".fts", ".FIT", ".FITS", ".FTS"}

METHOD_LABELS = {
    "skimage": "skimage phase correlation",
    "own": "own FFT phase correlation",
    "aa": "astroalign (translation only)",
    "aa_true": "astroalign (similarity)",
    "flow": "optical flow TV-L1",
}


@dataclass
class MethodResult:
    name: str
    dy: float = float("nan")
    dx: float = float("nan")
    rotation_deg: float = float("nan")
    scale: float = float("nan")
    aligned: np.ndarray | None = None
    residual: np.ndarray | None = None
    rms: float = float("nan")
    median_abs: float = float("nan")
    error: str | None = None
    flow_v: np.ndarray | None = field(default=None, repr=False)
    flow_u: np.ndarray | None = field(default=None, repr=False)


def _fits_in_dir(directory: Path) -> list[Path]:
    files = [
        p
        for p in sorted(directory.iterdir())
        if p.is_file() and p.suffix in FITS_SUFFIXES
    ]
    if not files:
        files = [
            p
            for p in sorted(directory.rglob("*"))
            if p.is_file() and p.suffix in FITS_SUFFIXES
        ]
    return files


def resolve_fits_path(path: str | Path, index: int) -> Path:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Not found: {p}")
    if p.is_file():
        return p
    files = _fits_in_dir(p)
    if not files:
        raise FileNotFoundError(f"No FITS files in {p}")
    if index < 0 or index >= len(files):
        raise IndexError(f"file index {index} out of range (n={len(files)}) in {p}")
    return files[index]


def _as_2d(data: np.ndarray) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    while arr.ndim > 2:
        arr = arr[0]
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2-D image, got shape {arr.shape}")
    return arr


def load_ccd(path: Path) -> CCDData:
    ccd = CCDData.read(path)
    ccd.data = _as_2d(ccd.data)
    if ccd.mask is None:
        ccd.mask = np.zeros(ccd.data.shape, dtype=bool)
    else:
        ccd.mask = np.asarray(ccd.mask, dtype=bool)
        if ccd.mask.shape != ccd.data.shape:
            ccd.mask = np.zeros(ccd.data.shape, dtype=bool)
    if ccd.uncertainty is None:
        ccd.uncertainty = StdDevUncertainty(np.ones_like(ccd.data, dtype=float))
    if getattr(ccd, "unit", None) is None:
        import astropy.units as u

        ccd.unit = u.adu
    return ccd


def downsample_ccd(ccd: CCDData, factor: int) -> CCDData:
    if factor <= 1:
        return ccd
    sl = (slice(None, None, factor), slice(None, None, factor))
    data = np.asarray(ccd.data, dtype=float)[sl]
    mask = np.asarray(ccd.mask, dtype=bool)[sl]
    unc = np.asarray(ccd.uncertainty.array, dtype=float)[sl]
    return CCDData(
        data,
        mask=mask,
        uncertainty=StdDevUncertainty(unc),
        meta=ccd.meta,
        unit=ccd.unit,
    )


def match_shapes(a: CCDData, b: CCDData) -> tuple[CCDData, CCDData]:
    ny = min(a.data.shape[0], b.data.shape[0])
    nx = min(a.data.shape[1], b.data.shape[1])
    if (a.data.shape[0], a.data.shape[1]) != (ny, nx):
        a = CCDData(
            a.data[:ny, :nx],
            mask=a.mask[:ny, :nx],
            uncertainty=StdDevUncertainty(a.uncertainty.array[:ny, :nx]),
            meta=a.meta,
            unit=a.unit,
        )
    if (b.data.shape[0], b.data.shape[1]) != (ny, nx):
        b = CCDData(
            b.data[:ny, :nx],
            mask=b.mask[:ny, :nx],
            uncertainty=StdDevUncertainty(b.uncertainty.array[:ny, :nx]),
            meta=b.meta,
            unit=b.unit,
        )
    return a, b


def _normalize(img: np.ndarray) -> np.ndarray:
    finite = img[np.isfinite(img)]
    if finite.size == 0:
        return img
    med = float(np.median(np.abs(finite)))
    if med <= 0:
        med = float(np.std(finite)) or 1.0
    return img / med


def residual_stats(reference: np.ndarray, aligned: np.ndarray) -> tuple[np.ndarray, float, float]:
    r = _normalize(reference)
    a = _normalize(aligned)
    res = a - r
    ok = np.isfinite(res)
    if not np.any(ok):
        return res, float("nan"), float("nan")
    return res, float(np.sqrt(np.mean(res[ok] ** 2))), float(np.median(np.abs(res[ok])))


def _write_temp_fits(ccd: CCDData, path: Path) -> None:
    ccd.write(path, overwrite=True)


def run_xy_method(name: str, ref: CCDData, other: CCDData, tmp: Path) -> MethodResult:
    from ost_photometry.reduce.registration import calculate_xy_image_shifts_core

    ref_path = tmp / "reference.fits"
    oth_path = tmp / "other.fits"
    _write_temp_fits(ref, ref_path)
    _write_temp_fits(other, oth_path)
    _id, shift, _flip = calculate_xy_image_shifts_core(
        str(oth_path), str(ref_path), 0, correlation_method=name
    )
    dy, dx = float(shift[0]), float(shift[1])
    result = MethodResult(name=name, dy=dy, dx=dx)
    if not np.isfinite(dy) or not np.isfinite(dx):
        result.error = "shift determination failed"
        return result
    aligned = ndi_shift(np.asarray(other.data, dtype=float), shift=(dy, dx), order=1)
    result.aligned = aligned
    result.residual, result.rms, result.median_abs = residual_stats(
        np.asarray(ref.data, dtype=float), aligned
    )
    return result


def run_aa_true(ref: CCDData, other: CCDData) -> MethodResult:
    from ost_photometry.reduce.registration import astro_align

    result = MethodResult(name="aa_true")
    aligned_ccd, transform = astro_align(ref, other)
    result.aligned = np.asarray(aligned_ccd.data, dtype=float)
    # skimage SimilarityTransform: translation is (x, y)
    result.dx = float(transform.translation[0])
    result.dy = float(transform.translation[1])
    result.rotation_deg = float(np.degrees(transform.rotation))
    result.scale = float(transform.scale)
    result.residual, result.rms, result.median_abs = residual_stats(
        np.asarray(ref.data, dtype=float), result.aligned
    )
    return result


def run_flow(ref: CCDData, other: CCDData) -> MethodResult:
    from skimage.registration import optical_flow_tvl1

    from ost_photometry.reduce.registration import optical_flow_align

    result = MethodResult(name="flow")
    aligned_ccd = optical_flow_align(ref, other)
    result.aligned = np.asarray(aligned_ccd.data, dtype=float)
    flow_v, flow_u = optical_flow_tvl1(
        np.asarray(ref.data, dtype=float),
        np.asarray(other.data, dtype=float),
    )
    result.flow_v = flow_v
    result.flow_u = flow_u
    result.dy = float(np.nanmedian(flow_v))
    result.dx = float(np.nanmedian(flow_u))
    result.residual, result.rms, result.median_abs = residual_stats(
        np.asarray(ref.data, dtype=float), result.aligned
    )
    return result


def run_method(name: str, ref: CCDData, other: CCDData, tmp: Path) -> MethodResult:
    try:
        if name in {"skimage", "own", "aa"}:
            return run_xy_method(name, ref, other, tmp)
        if name == "aa_true":
            return run_aa_true(ref, other)
        if name == "flow":
            return run_flow(ref, other)
        raise ValueError(f"Unknown method {name!r}")
    except Exception as exc:
        return MethodResult(name=name, error=f"{type(exc).__name__}: {exc}")


def _zscale_imshow(ax, data: np.ndarray, *, cmap: str = "gray") -> None:
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        ax.set_axis_off()
        return
    try:
        norm = ImageNormalize(data, interval=ZScaleInterval())
    except Exception:
        lo, hi = np.percentile(finite, (1.0, 99.5))
        norm = ImageNormalize(vmin=lo, vmax=max(hi, lo + 1e-6))
    ax.imshow(data, origin="lower", cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_shift_summary(results: list[MethodResult], path: Path, *, scale: float) -> None:
    names = [r.name for r in results]
    dx = np.array([r.dx * scale for r in results])
    dy = np.array([r.dy * scale for r in results])
    rms = np.array([r.rms for r in results])
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    x = np.arange(len(names))
    axes[0].bar(x - 0.18, dx, width=0.36, label="Δx")
    axes[0].bar(x + 0.18, dy, width=0.36, label="Δy")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(names, rotation=30, ha="right")
    axes[0].set_ylabel("Shift [original pix]")
    axes[0].legend(fontsize=8)
    axes[0].set_title("Translation")
    axes[0].grid(True, axis="y", alpha=0.3)

    rot = np.array([r.rotation_deg for r in results])
    scl = np.array([r.scale for r in results])
    axes[1].scatter(x, rot, c="C2", label="rotation [deg]")
    ax1b = axes[1].twinx()
    ax1b.scatter(x, scl, c="C3", marker="s", label="scale")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(names, rotation=30, ha="right")
    axes[1].set_ylabel("Rotation [deg]")
    ax1b.set_ylabel("Scale")
    axes[1].set_title("Similarity (aa_true)")
    axes[1].grid(True, axis="y", alpha=0.3)

    axes[2].bar(x, rms, color="0.45")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(names, rotation=30, ha="right")
    axes[2].set_ylabel("RMS residual (norm. ADU)")
    axes[2].set_title("Alignment residual")
    axes[2].grid(True, axis="y", alpha=0.3)
    fig.suptitle("Registration method comparison", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_shift_differences(results: list[MethodResult], path: Path, *, scale: float) -> None:
    n = len(results)
    names = [r.name for r in results]
    dxy = np.full((n, n), np.nan)
    for i, a in enumerate(results):
        for j, b in enumerate(results):
            if not (np.isfinite(a.dx) and np.isfinite(b.dx)):
                continue
            dxy[i, j] = scale * np.hypot(a.dx - b.dx, a.dy - b.dy)
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    im = ax.imshow(dxy, origin="upper", cmap="magma")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=30, ha="right")
    ax.set_yticklabels(names)
    for i in range(n):
        for j in range(n):
            val = dxy[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="w")
    fig.colorbar(im, ax=ax, label="|Δshift| [original pix]")
    ax.set_title("Pairwise translation difference")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_method_panels(
    ref: np.ndarray,
    other: np.ndarray,
    result: MethodResult,
    path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 8.6))
    _zscale_imshow(axes[0, 0], ref)
    axes[0, 0].set_title("Reference")
    _zscale_imshow(axes[0, 1], other)
    axes[0, 1].set_title("Other (unaligned)")
    if result.aligned is not None:
        _zscale_imshow(axes[1, 0], result.aligned)
    axes[1, 0].set_title("Other aligned")
    if result.residual is not None:
        vmax = np.nanpercentile(np.abs(result.residual), 98)
        axes[1, 1].imshow(
            result.residual,
            origin="lower",
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest",
        )
        axes[1, 1].set_xticks([])
        axes[1, 1].set_yticks([])
    axes[1, 1].set_title("Aligned − reference (norm.)")
    extra = ""
    if np.isfinite(result.rotation_deg):
        extra = f", rot={result.rotation_deg:.3f}°, scale={result.scale:.5f}"
    err = f"\nFAILED: {result.error}" if result.error else ""
    fig.suptitle(
        f"{METHOD_LABELS.get(result.name, result.name)}\n"
        f"Δx={result.dx:.3f} pix, Δy={result.dy:.3f} pix, "
        f"RMS={result.rms:.4f}{extra}{err}",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_residual_grid(results: list[MethodResult], path: Path) -> None:
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 3.6), squeeze=False)
    for ax, res in zip(axes[0], results):
        if res.residual is None:
            ax.set_axis_off()
            ax.set_title(res.name)
            continue
        vmax = np.nanpercentile(np.abs(res.residual), 98)
        ax.imshow(
            res.residual,
            origin="lower",
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{res.name}\nRMS={res.rms:.3f}")
    fig.suptitle("Residuals after alignment (same stretch per panel)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_flow_field(result: MethodResult, path: Path) -> None:
    if result.flow_u is None or result.flow_v is None:
        return
    ny, nx = result.flow_u.shape
    step = max(1, min(ny, nx) // 25)
    yy, xx = np.mgrid[0:ny:step, 0:nx:step]
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    mag = np.hypot(result.flow_u, result.flow_v)
    im = ax.imshow(mag, origin="lower", cmap="viridis")
    ax.quiver(
        xx,
        yy,
        result.flow_u[::step, ::step],
        result.flow_v[::step, ::step],
        color="w",
        angles="xy",
        scale_units="xy",
        scale=1.0,
        width=0.002,
        alpha=0.8,
    )
    fig.colorbar(im, ax=ax, label="|flow| [pix]")
    ax.set_title(
        f"Optical flow (median Δx={result.dx:.3f}, Δy={result.dy:.3f} pix)"
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_csv(results: list[MethodResult], path: Path, *, scale: float) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "method",
                "dx_original_pix",
                "dy_original_pix",
                "rotation_deg",
                "scale",
                "rms_residual",
                "median_abs_residual",
                "error",
            ]
        )
        for r in results:
            writer.writerow(
                [
                    r.name,
                    r.dx * scale if np.isfinite(r.dx) else "",
                    r.dy * scale if np.isfinite(r.dy) else "",
                    r.rotation_deg if np.isfinite(r.rotation_deg) else "",
                    r.scale if np.isfinite(r.scale) else "",
                    r.rms if np.isfinite(r.rms) else "",
                    r.median_abs if np.isfinite(r.median_abs) else "",
                    r.error or "",
                ]
            )


def reduce_dataset(dataset: Path, output: Path) -> Path:
    from ost_photometry.calibration_parameters import get_image_types
    from ost_photometry.reduce.redu import reduce_main

    output.mkdir(parents=True, exist_ok=True)
    reduce_main(
        str(dataset),
        str(output),
        image_type_dir=get_image_types(),
        stack_images=False,
        find_wcs=False,
        save_only_transformation=True,
        estimate_fwhm=False,
        debug=False,
    )
    light = output / "light"
    if not light.is_dir():
        raise FileNotFoundError(f"reduce_main did not write {light}")
    return light


def synthetic_pair(shift_yx: tuple[float, float] = (2.5, -3.0)) -> tuple[CCDData, CCDData]:
    rng = np.random.default_rng(0)
    ny, nx = 180, 220
    yy, xx = np.mgrid[0:ny, 0:nx]
    ref = np.zeros((ny, nx), dtype=float)
    for _ in range(18):
        y0 = rng.uniform(15, ny - 15)
        x0 = rng.uniform(15, nx - 15)
        amp = rng.uniform(80.0, 250.0)
        ref += amp * np.exp(-((yy - y0) ** 2 + (xx - x0) ** 2) / (2.3**2))
    ref += 20.0 + rng.normal(0.0, 1.5, size=ref.shape)
    other = ndi_shift(ref, shift=shift_yx, order=3)
    try:
        import astropy.units as u

        unit = u.adu
    except Exception:
        unit = "adu"
    kw = dict(
        mask=np.zeros_like(ref, dtype=bool),
        uncertainty=StdDevUncertainty(np.ones_like(ref)),
        unit=unit,
    )
    return CCDData(ref, **kw), CCDData(other, **kw)


def compare_frames(
    ref: CCDData,
    other: CCDData,
    output: Path,
    *,
    methods: list[str],
    downsample: int,
    plot_format: str,
) -> list[MethodResult]:
    output.mkdir(parents=True, exist_ok=True)
    ref, other = match_shapes(ref, other)
    ref = downsample_ccd(ref, downsample)
    other = downsample_ccd(other, downsample)
    scale = float(max(downsample, 1))
    results: list[MethodResult] = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for name in methods:
            print(f"  {name} …", flush=True)
            results.append(run_method(name, ref, other, tmp_path))
            if results[-1].error:
                print(f"    failed: {results[-1].error}", flush=True)
            else:
                print(
                    f"    Δx={results[-1].dx * scale:.3f}  "
                    f"Δy={results[-1].dy * scale:.3f} pix (original)",
                    flush=True,
                )

    ref_data = np.asarray(ref.data, dtype=float)
    other_data = np.asarray(other.data, dtype=float)
    plot_shift_summary(results, output / f"registration_shift_summary.{plot_format}", scale=scale)
    plot_shift_differences(
        results, output / f"registration_shift_differences.{plot_format}", scale=scale
    )
    plot_residual_grid(results, output / f"registration_residuals.{plot_format}")
    for res in results:
        plot_method_panels(
            ref_data,
            other_data,
            res,
            output / f"registration_{res.name}.{plot_format}",
        )
        if res.name == "flow":
            plot_flow_field(res, output / f"registration_flow_field.{plot_format}")
    write_csv(results, output / "registration_shifts.csv", scale=scale)
    return results


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare ost_photometry registration methods on two frames."
    )
    p.add_argument("--ref", default=PATH_REFERENCE, help="Reference FITS or directory")
    p.add_argument("--other", default=PATH_OTHER, help="Second FITS or directory")
    p.add_argument("--out", default=PATH_OUTPUT, help="Output directory")
    p.add_argument("--ref-index", type=int, default=FILE_INDEX_REFERENCE)
    p.add_argument("--other-index", type=int, default=FILE_INDEX_OTHER)
    p.add_argument("--reduce", action="store_true", default=REDUCE_FIRST)
    p.add_argument("--downsample", type=int, default=DOWNSAMPLE)
    p.add_argument(
        "--methods",
        nargs="+",
        default=METHODS,
        choices=list(METHOD_LABELS),
    )
    p.add_argument("--format", dest="plot_format", default=PLOT_FORMAT)
    p.add_argument(
        "--self-test",
        action="store_true",
        help="Run on a synthetic shifted star field instead of FITS files",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out = Path(args.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    if args.self_test:
        print("Self-test: synthetic field, true shift Δx=-3.0, Δy=+2.5 pix")
        ref, other = synthetic_pair((2.5, -3.0))
        compare_frames(
            ref,
            other,
            out,
            methods=list(args.methods),
            downsample=max(1, args.downsample),
            plot_format=args.plot_format,
        )
        print(f"Wrote plots to {out}")
        return 0

    if args.ref in {"?", ""} or args.other in {"?", ""}:
        print("Set PATH_REFERENCE / PATH_OTHER in the script, or pass --ref and --other.")
        return 2

    ref_path = Path(args.ref).expanduser()
    other_path = Path(args.other).expanduser()
    if args.reduce:
        print(f"Reducing {ref_path} …")
        ref_dir = reduce_dataset(ref_path, out / "reduce_reference")
        print(f"Reducing {other_path} …")
        other_dir = reduce_dataset(other_path, out / "reduce_other")
        ref_fits = resolve_fits_path(ref_dir, args.ref_index)
        other_fits = resolve_fits_path(other_dir, args.other_index)
    else:
        ref_fits = resolve_fits_path(ref_path, args.ref_index)
        other_fits = resolve_fits_path(other_path, args.other_index)

    print(f"Reference: {ref_fits}")
    print(f"Other:     {other_fits}")
    ref = load_ccd(ref_fits)
    other = load_ccd(other_fits)
    compare_frames(
        ref,
        other,
        out,
        methods=list(args.methods),
        downsample=max(1, args.downsample),
        plot_format=args.plot_format,
    )
    print(f"Wrote plots to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
