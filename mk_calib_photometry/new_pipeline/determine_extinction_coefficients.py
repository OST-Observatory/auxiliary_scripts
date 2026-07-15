#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
PURPOSE:
     Determine atmospheric extinction coefficients from reduced observations
     using the cat-star.org method (flux/magnitude vs airmass).

EXPLANATION:
     Observe stars (e.g. G2V) over several hours as they rise or set.
     The pipeline runs WCS, extraction, and correlation, then fits
     extinction coefficients k per filter. Output: extinction_coefficients.json
     and context.extinction_coefficients for use in differential calibration.

     After several nights, aggregate JSON files into the site table:
       python scripts/aggregate_site_extinction.py --nights ... --out ...
     See ost_photometry_package/docs/EXTINCTION_COEFFICIENTS.md

     Requires multiple images per filter at different airmasses.
     Best done on clear nights; moonlit nights are acceptable.
"""

############################################################################
#             Configuration: modify the file in this section               #
############################################################################

###
#   Cluster or field identifier (e.g., NGC 381, or "extinction_field")
#
name: str = "?"

###
#   Finder options
#
fwhm: float | None = None

###
#   Aperture or ePSF photometry
#
photometry_extraction_method: str = "APER"
# photometry_extraction_method: str = "PSF"

###
#   Magnitude column for extinction fit (default: mags_fit from extraction)
#
extinction_fit_mag_col: str = "mags_fit"

#   Use flux instead of magnitude (e.g. flux_fit)
# extinction_fit_use_flux: bool = False

############################################################################
#                Additional options: only edit if necessary                #
#              Add up to 10 variable set for different filter              #
############################################################################

###
# Filter 1 (e.g., B, V, R,...)
#
filter_1: str = "B"

# Path to the images of filter 1 (directory with multiple images, or single file)
# For extinction fit: multiple images at different airmasses required!
path_1: str = "output/filter_B/"

###
# Filter 2 (e.g., B, V, R,...)
#
filter_2: str = "V"

# Path to the images of filter 2
path_2: str = "output/filter_V/"

# Path to store the output
output_dir: str = "output/"

# Output filename for extinction coefficients
extinction_coefficients_filename: str = "extinction_coefficients.json"

############################################################################
#                               Libraries                                  #
############################################################################

import json
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from ost_photometry import style
from ost_photometry.analyze import analyze
from ost_photometry.analyze.pipeline import PipelineConfig

############################################################################
#                                  Main                                    #
############################################################################

if __name__ == "__main__":
    #   Set start time
    start_time = time.time()

    #   Prepare variable lists and dictionaries from the individual
    #   definitions above
    filter_list: list[str] = []
    img_paths: dict[str, str] = {}
    fwhm_object_psf: dict[str, float] = {}
    for i in range(0, 10):
        if "filter_" + str(i) in locals():
            filter_list.append(locals()["filter_" + str(i)])
            img_paths[locals()["filter_" + str(i)]] = locals()["path_" + str(i)]
            fwhm_object_psf[locals()["filter_" + str(i)]] = fwhm

    ###
    #   Pipeline config: run WCS, extraction, correlation, extinction fit;
    #   skip calibration (post-process steps skip when calibration skipped)
    #
    config = PipelineConfig(
        extinction_mode="from_value_airmass",
        skip_calibration=True,
        protect_calibration_objects=True,
        photometry_extraction_method=photometry_extraction_method,
        extinction_fit_mag_col=extinction_fit_mag_col,
        extinction_coefficients_filename=extinction_coefficients_filename,
        fwhm_object_psf=fwhm_object_psf,
    )

    ###
    #   Initialize observation and run pipeline
    #
    observation = analyze.Observation(object_names=[name])
    observation.run_pipeline(
        filter_list,
        image_paths=img_paths,
        output_dir=output_dir,
        config=config,
        extraction_mode="auto",
    )

    ###
    #   Report result (load from saved JSON)
    #
    out_path = Path(output_dir) / extinction_coefficients_filename
    if out_path.exists():
        with open(out_path) as f:
            coeffs = json.load(f)
        print(style.Bcolors.OKGREEN + "   Extinction coefficients:" + style.Bcolors.ENDC)
        for filt, c in coeffs.items():
            k = c.get("k_prime", "?")
            e = c.get("k_prime_err", "?")
            print(f"     {filt}: k' = {k} ± {e} mag/airmass")

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))

"""
    Change Log
    ----------
        0.1   (2025)
           - initial release
           - extinction fit via cat-star.org method (flux/mag vs airmass)
"""
