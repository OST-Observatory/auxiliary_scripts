#! /usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
####          Configuration: modify the file in this section            ####
############################################################################

###
#   Cluster identifier
#
nameobj = '?'


########################  Filter 1  #########################
### Define filter (e.g., U, B,V,...)
filter_1 = '?'

### Path to the images of filter 1
img_1    = '?'


########################  Filter 2  #########################
### Define filter 2 (e.g., U, B,V,...)
filter_2 = '?'

### Path to the images of filter 1
img_2    = '?'

###
#   Finder options
#
# Set sigma -> characterizes the size of the diffraction patterns
sigma = 3.0


############################################################################
####             Additional options: only edit if necessary             ####
############################################################################


###
#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
#
outdir='output/'

###
#   Aperture or ePSF photometry
#
#   APER or PSF
photometry = 'APER'
#photometry = 'PSF'


###
#   Valid filter combinations to calculate magnitude transformation
#   dict -> key = filter, value = list(first color, second color)
#
valid_calibs = [['U','V'], ['B','V'], ['V','R'], ['V','I']]


###
#   Calibration source (possibilities: simbad_vot, UCAC4, GSC2.3, URAT1,
#                                      NOMAD, APASS)
#
#calib_methode = 'simbad_vot' # Downloaded table in VO format from Simbad
calib_methode = 'vsp'        # Load calibration data for variable stars
                             # directly from the AAVSO website
#calib_methode = 'UCAC4'      # UCAC4 Catalogue (Zacharias+, 2012)
#calib_methode = 'GSC2.3'     # The Full GSC2.3.2 Catalogue
#calib_methode = 'URAT1'      # URAT1 Catalog (Zacharias+ 2015)
#calib_methode = 'NOMAD'      # NOMAD Catalog (Zacharias+ 2005)
calib_methode = 'APASS'      # AAVSO Photo. All Sky Surv. DR9(Henden+,2016)

#   Dictionary with catalog information
vizier_dict = {
    'UCAC4':'I/322A',
    'GSC2.3':'I/305',
    'URAT1':'I/329',
    'NOMAD':'I/297',
    'APASS':'II/336/apass9',
}

#   File with calibration stars (has to be in VO format)
file_calib = None

#   Magnitude range of the calibration stars
mag_range = (12., 15.)

#   Apply weights in calculation of magnitude transformation
weights = True


###
#   Additional finder options
#
#   Extraction methode: DAO or IRAF
methode = 'DAO'
#methode = 'IRAF'

#   Sigma background
sigma_bkg = 3

## Threshold multiplier:
#   First iteration:
multi_start = 7.

#   Final iteration:
multi = 3.
multi = 5.

#   DAO grouper options
multi_grouper = 2.0

#   Remove objects with negative flux uncertainties
strict_cleaning = True
#strict_cleaning = False


###
#   ePSF options
#
#   Oversampling
oversampling = 2
#oversampling = 4

#   Max. number of iterations
maxiters = 12
maxiters = 7
#maxiters = 3

#   Size extraction box (pixel)
size_epsf = 25

#   Fraction of all stars used for EPSF determination
frac_epsf_stars = 0.05
frac_epsf_stars = 0.1
frac_epsf_stars = 0.15
frac_epsf_stars = 0.2

#   Minimal number of required ePSF stars
min_eps_stars = 25

#   Require that 'min_eps_stars' will be reached
#strict_eps=False
strict_eps=True


###
#   Aperture options
#
#   Extraction radius stars in arcsec or pixel
rstars = 4.
#rstars = 5.
#   Extraction radius background (inner and outer radii) in arcsec or pixel
rbg_in  = 7.
rbg_out = 10.
#   Unit
#r_unit = 'pixel'
r_unit = 'arcsec'


###
#   newsrcor options
#
#   ID of the reference image
ref_ID = 0

#   Critical radius outside which correlations are rejected
dcr     = 3
#dcr     = 5
#dcr     = 7

#   Refinement option (3: take all stars within the distance of 'dcr' to
#   another as a match, 0: take only the one with the minimal distance, 1:
#   force in addition the identification to be one-to-one, 2: more rigorous
#   version of 1)
option = 1

#   Fraction of bad images - Used to reject bad objects:
#   Objects that are not on the 'bfrac' fraction of all images will
#   be rejected (works counter intuitively, because badd sources often lead
#   to the rejection of good images => small bfrac values might result in a
#   lot of recjected images)
bfrac   = 0.8
bfrac   = 0.9


###
#   WCS options
#
#   Methode to determine WCS
wcs_method = 'astrometry'                  #   -> astrometry.net
#wcs_method = 'twirl'                       #   -> twirl libary


###
#   Plot options
#
#   Make star map for the initial extraction
plot_ifi = True

#   Make only the star map plot for the reference image [refid]
plot_test = True


###
#   Cosmic ray removal
#   Note: For cosmic ray removal, images should be pre-reduced with
#   1_reduce_images.py before running this script.
#
rmcos = True
rmcos = False

#   Parameters (used when running with pre-reduced images):
objlim  = 5.
sigclip = 4.0

#   Needs to changed in the future
camera = 'QHY600M'


###
#   Multiprocessing
#
ncores = 6


###
#   Expert interface
#   -> add here, if images from more than 2 filters should be reduced
#
filter_list = [filter_1, filter_2]
img_dirs    = {filter_1:img_1, filter_2:img_2}
sigma_psf   = {filter_1:sigma, filter_2:sigma}

###
#   Camera specific parameters (for reference; cosmic ray removal
#   is typically done in 1_reduce_images.py)
#
if camera == 'STF8300':
    readnoise = 9.3
    gain      = None
    dark_rate = {0:0.18, -10:0.04, -15.8:0.02}
    satlevel  = 65535.
elif camera == 'QHY600M':
    readnoise = 7.904
    gain      = 1.292
    dark_rate = {-20:0.0022, -10:0.0046}
    satlevel  = 65535.
else:
    raise RuntimeError(
        "Error: camera type not known\n"
        "\t-> check variable: camera\n"
        "\t-> Exit\n"
    )


############################################################################
####                            Libraries                               ####
############################################################################

from os.path import join

import time

import warnings
warnings.filterwarnings('ignore')

from astropy.table import Table
from ost_photometry import checks
from ost_photometry import style
from ost_photometry.analyze import Observation
from ost_photometry.analyze.models import ImageSeries
from ost_photometry.analyze import calibration, correlate, utilities
from ost_photometry.analyze.extraction import extract_multiprocessing


############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Set start time
    start_time = time.time()

    ###
    #   Initialize observation container
    #
    observation = Observation()

    ###
    #   Check output directories
    #
    checks.check_output_directories(
        outdir,
        join(outdir, 'tables'),
    )

    ###
    #   Check image directories
    #
    checks.check_dir(img_dirs)

    #   Outer loop over all filter
    for filt in filter_list:
        print(
            style.Bcolors.HEADER
            + "   Analyzing " + filt + " images"
            + style.Bcolors.ENDC
        )

        #   Initialize image series object
        observation.image_series_dict[filt] = ImageSeries(
            filt,
            img_dirs[filt],
            outdir,
            ref_ID,
        )

        ###
        #   Find the WCS solution for the image
        #
        utilities.find_wcs(
            observation.image_series_dict[filt],
            reference_image_index=ref_ID,
            method=wcs_method,
            indent=2,
        )

        ###
        #   Main extraction of object positions and object fluxes
        #   using multiprocessing
        #
        extract_multiprocessing(
            observation.image_series_dict[filt],
            ncores,
            fwhm_object_psf=sigma_psf,
            sigma_value_background_clipping=sigma_bkg,
            multiplier_background_rms=multi_start,
            size_epsf_region=size_epsf,
            fraction_epsf_stars=frac_epsf_stars,
            oversampling_factor_epsf=oversampling,
            max_n_iterations_epsf_determination=maxiters,
            object_finder_method=methode,
            multiplier_background_rms_epsf=multi,
            multiplier_grouper_epsf=multi_grouper,
            strict_cleaning_epsf_results=strict_cleaning,
            minimum_n_eps_stars=min_eps_stars,
            strict_epsf_checks=strict_eps,
            photometry_extraction_method=photometry,
            radius_aperture=rstars,
            inner_annulus_radius=rbg_in,
            outer_annulus_radius=rbg_out,
            radii_unit=r_unit,
            plots_for_all_images=plot_ifi,
        )

        ###
        #   Correlate results from all images, while preserving the
        #   calibration stars
        #
        correlate.correlate_preserve_calibration_objects(
            observation.image_series_dict[filt],
            filter_list,
            calibration_source=calib_methode,
            calibration_catalog_mag_range=mag_range,
            vizier_dict=vizier_dict,
            calib_file=file_calib,
            max_pixel_between_objects=dcr,
            ooi_correlation_strategy=option,
            cross_identification_limit=1,
            reference_image_index=ref_ID,
            n_allowed_non_detections_object=1,
            expected_bad_image_fraction=bfrac,
            protect_calibration_objects=False,
            plot_only_reference_starmap=plot_test,
        )

    ###
    #   Make new calibration table and add object name to
    #   the calibration table
    #
    tbl_trans            = Table()
    tbl_trans['name']    = [nameobj]

    #   Loop over allowed filter combinations to allow for the calculation
    #   of the transformation coefficients
    for calib_fil in valid_calibs:
        #   Check if filter combination is valid
        if calib_fil[0] in filter_list and calib_fil[1] in filter_list:
            for i in range(0, len(calib_fil)):
                key = calib_fil[i]

                #   Set up filter list
                filt_list    = [calib_fil[0], calib_fil[1]]

                #   Add air mass and object to the calibration table
                tbl_trans['airmass_'+key] = [
                    observation.image_series_dict[key].median_air_mass()
                ]

                ###
                #   Correlate the results from the different filter and
                #   determine transformation coefficients
                #
                calibration.calculate_trans(
                    observation,
                    key,
                    filt_list,
                    tbl_trans,
                    apply_uncertainty_weights=weights,
                    max_pixel_between_objects=dcr,
                    ooi_correlation_strategy=option,
                    calibration_source=calib_methode,
                    vizier_dict=vizier_dict,
                    calibration_file=file_calib,
                    calibration_catalog_mag_range=mag_range,
                )

                tbl_trans['jd'] = [
                    observation.image_series_dict[key].median_observation_time()
                ]

    #   Write table and check output directories
    tbl_trans.write(
        outdir+'/tables/trans_para_'+nameobj.replace(' ', '_')+'.dat',
        format='ascii',
        overwrite=True,
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
