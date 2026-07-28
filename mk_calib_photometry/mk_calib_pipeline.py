"""Shared helpers for mk_calib_photometry scripts using ``Observation.run_pipeline``."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ost_photometry.analyze.calibration.mk_calib import (
    calibrate_mk_calib_filter_pair,
    merge_field_transformation_records,
    write_field_transformation_json,
    write_trans_para_table,
)
from ost_photometry.analyze.pipeline import PipelineConfig

if TYPE_CHECKING:
    from ost_photometry.analyze import Observation


def build_transformation_pipeline_config(
    *,
    photometry: str,
    wcs_method: str,
    ref_id: int,
    sigma_psf: dict[str, float],
    ncores: int,
    methode: str,
    sigma_bkg: float,
    multi_start: float,
    multi: float,
    multi_grouper: float,
    strict_cleaning: bool,
    oversampling: int,
    maxiters: int,
    size_epsf: int,
    frac_epsf_stars: float,
    min_eps_stars: int,
    strict_eps: bool,
    rstars: float,
    rbg_in: float,
    rbg_out: float,
    r_unit: str,
    plot_ifi: bool,
    dcr: int,
    option: int,
    maxid: int,
    nmissed: int,
    bfrac: float,
    protect_calibration_objects: bool,
    calib_methode: str,
    vizier_dict: dict[str, str],
    file_calib: str | None,
    mag_range: tuple[float, float],
    rmcos: bool = False,
    readnoise: float = 8.0,
    sigclip: float = 4.5,
    satlevel: float = 65535.0,
    objlim: float = 5.0,
) -> PipelineConfig:
    """Pipeline config for WCS, extraction, and intra correlation (mk_calib step 2)."""
    return PipelineConfig.from_preset(
        "extract_protect_calibrators",
        overrides={
            "wcs_method": wcs_method,
            "photometry_extraction_method": photometry,
            "reference_image_index": ref_id,
            "fwhm_object_psf": sigma_psf,
            "n_cores_multiprocessing": ncores,
            "object_finder_method": methode,
            "sigma_value_background_clipping": sigma_bkg,
            "multiplier_background_rms": multi_start,
            "multiplier_background_rms_epsf": multi,
            "multiplier_grouper_epsf": multi_grouper,
            "strict_cleaning_epsf_results": strict_cleaning,
            "oversampling_factor_epsf": oversampling,
            "max_n_iterations_epsf_determination": maxiters,
            "size_epsf_region": size_epsf,
            "fraction_epsf_stars": frac_epsf_stars,
            "minimum_n_eps_stars": min_eps_stars,
            "strict_epsf_checks": strict_eps,
            "radius_aperture": rstars,
            "inner_annulus_radius": rbg_in,
            "outer_annulus_radius": rbg_out,
            "radii_unit": r_unit,
            "plots_for_all_images": plot_ifi,
            "max_pixel_between_objects": dcr,
            "ooi_correlation_strategy": option,
            "cross_identification_limit": maxid,
            "n_allowed_non_detections_object": nmissed,
            "expected_bad_image_fraction": bfrac,
            "protect_calibration_objects": protect_calibration_objects,
            "calibration_source": calib_methode,
            "vizier_dict": vizier_dict,
            "path_calibration_file": file_calib,
            "calibration_catalog_mag_range": mag_range,
            "cosmic_ray_removal": rmcos,
            "read_noise": readnoise,
            "sigma_clipping_value": sigclip,
            "saturation_level": satlevel,
            "limiting_contrast_rm_cosmics": objlim,
        },
    )


def build_calibration_pipeline_config(
    extraction_config: PipelineConfig,
    *,
    apply_weights: bool = True,
) -> PipelineConfig:
    """Config for per-pair CalibrationEngine calibration after extraction."""
    flat = extraction_config.as_flat_dict()
    cal_keys = (
        "calibration_source",
        "vizier_dict",
        "path_calibration_file",
        "calibration_catalog_mag_range",
        "max_pixel_between_objects",
        "ooi_correlation_strategy",
        "cross_identification_limit",
        "n_allowed_non_detections_object",
        "expected_bad_image_fraction",
        "protect_calibration_objects",
        "protect_ooi",
        "protected_object_ids",
        "correlation_method",
        "separation_limit",
        "reference_image_index",
        "verbose",
        "file_type_plots",
        "distribution_samples",
        "zp_subsample_statistic",
    )
    overrides = {k: flat[k] for k in cal_keys if k in flat}
    overrides["skip_calibration"] = False
    return PipelineConfig.from_preset("linear_fit_ensemble", overrides=overrides)


def run_transformation_extraction_pipeline(
    observation: "Observation",
    filter_list: list[str],
    img_dirs: dict[str, str],
    outdir: str,
    config: PipelineConfig,
) -> None:
    """Run WCS, extraction, and intra-filter correlation via the standard pipeline."""
    observation.run_pipeline(
        filter_list,
        image_paths=img_dirs,
        output_dir=outdir,
        config=config,
        extraction_mode="multi",
    )


def write_field_transformation_table(
    observation: "Observation",
    *,
    nameobj: str,
    filter_list: list[str],
    valid_calibs: list[list[str]],
    outdir: str,
    extraction_config: PipelineConfig,
    weights: bool = True,
) -> tuple[str, str]:
    """
    Calibrate each filter pair with CalibrationEngine and write field tables.

    Writes legacy ``trans_para_*.dat`` and JSON sidecar for second-order analysis.

    Returns ``(ascii_path, json_path)``.
    """
    cal_config = build_calibration_pipeline_config(extraction_config)
    records = []
    for calib_fil in valid_calibs:
        if calib_fil[0] not in filter_list or calib_fil[1] not in filter_list:
            continue
        pair = [calib_fil[0], calib_fil[1]]
        records.append(
            calibrate_mk_calib_filter_pair(
                observation,
                pair,
                cal_config,
                apply_weights=weights,
            )
        )
    if not records:
        raise RuntimeError("No valid filter pairs for mk_calib calibration")

    merged = merge_field_transformation_records(records)
    merged.name = nameobj

    stem = f"{outdir}/tables/trans_para_{nameobj.replace(' ', '_')}"
    ascii_path = write_trans_para_table(merged, f"{stem}.dat")
    json_path = write_field_transformation_json(merged, f"{stem}.json")
    return str(ascii_path), str(json_path)
