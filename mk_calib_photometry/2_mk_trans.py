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
filter_1   = '?'

### Path to the images of filter 1
dir_1      = '?'


########################  Filter 2  #########################
### Define filter 2 (e.g., U, B,V,...)
filter_2   = '?'

### Path to the images of filter 1
dir_2      = '?'


###
#   Finder options
#
# Set sigma -> characterizes the size of the diffraction patterns
sigma = 3.


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
#   aper or PSF
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
calib_methode = 'ASCII'      # Load ASCII table -> set col_dict dictonary
calib_methode = 'vsp'        # Load calibration data for variable stars
                             # directly from the AAVSO website
#calib_methode = 'UCAC4'      # UCAC4 Catalogue (Zacharias+, 2012)
#calib_methode = 'GSC2.3'     # The Full GSC2.3.2 Catalogue
calib_methode = 'URAT1'       # URAT1 Catalog (Zacharias+ 2015)
#calib_methode = 'NOMAD'      # NOMAD Catalog (Zacharias+ 2005)
#calib_methode = 'APASS'      # AAVSO Photo. All Sky Surv. DR9(Henden+,2016)

#   Dictionary with catalog information
vizier_dict = {'UCAC4':'I/322A', 'GSC2.3':'I/305', 'URAT1':'I/329',
               'NOMAD':'I/297', 'APASS':'II/336'}

#   File with calibration stars
file_calib = None

#   Magnitude limit of the calibration stars
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
multi = 7.

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

#   Max. number of iterations
maxiters = 12
#maxiters = 7
maxiters = 5
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
ref_ID   = 1
#ref_ID   = 2
#ref_ID   = 6

#   Critical radius outside which correlations are rejected
dcr     = 3
#dcr     = 5
dcr     = 7
dcr     = 9
#dcr     = 13
#dcr     = 20

#   Max. number of identical identifications, images with higher values will
#   be rejected
maxid   = 1

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
bfrac   = 0.4
bfrac   = 0.8
bfrac   = 0.9
#bfrac   = 1.

#   Limit for the number of images on which an object is not found.
#   When this limit is reached, the corresponding object is discarded.
nmissed = 1
nmissed = 3
#nmissed = 10

#   Preserve calibration stars in cross correlation between different images
s_refOBJ = False
#s_refOBJ = True


###
#   WCS options
#
#   Create a WCS solution for all images
mk_all_wcs = False
#   Methode to determine WCS
wcs_method = 'astrometry'                  #   -> astrometry.net
#wcs_method = 'twirl'                       #   -> twirl libary


###
#   Plot options
#
#   Make star map plots for all stars
plot_ifi = False
#plot_ifi = True

#   Make only the star map plot for the reference image [refid]
plot_test = True

###
#   Multiprocessing
#
ncores = 6


###
#   Verbose output
#
verbose = False


###
#   Expert interface
#   -> add here, if images from more than 2 filters should be reduced
#
filter_list = [filter_1, filter_2]
img_dirs    = {filter_1:dir_1, filter_2:dir_2}
sigma_psf   = {filter_1:sigma, filter_2:sigma}


############################################################################
####                            Libraries                               ####
############################################################################

import time

import warnings
warnings.filterwarnings('ignore')

from ost_photometry import style
from ost_photometry.analyze import Observation

from mk_calib_pipeline import (
    build_transformation_pipeline_config,
    run_transformation_extraction_pipeline,
    write_field_transformation_table,
)


############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    start_time = time.time()

    observation = Observation(object_names=[nameobj])

    pipeline_config = build_transformation_pipeline_config(
        photometry=photometry,
        wcs_method=wcs_method,
        ref_id=ref_ID,
        sigma_psf=sigma_psf,
        ncores=ncores,
        methode=methode,
        sigma_bkg=sigma_bkg,
        multi_start=multi_start,
        multi=multi,
        multi_grouper=multi_grouper,
        strict_cleaning=strict_cleaning,
        oversampling=oversampling,
        maxiters=maxiters,
        size_epsf=size_epsf,
        frac_epsf_stars=frac_epsf_stars,
        min_eps_stars=min_eps_stars,
        strict_eps=strict_eps,
        rstars=rstars,
        rbg_in=rbg_in,
        rbg_out=rbg_out,
        r_unit=r_unit,
        plot_ifi=plot_ifi,
        dcr=dcr,
        option=option,
        maxid=maxid,
        nmissed=nmissed,
        bfrac=bfrac,
        protect_calibration_objects=s_refOBJ,
        calib_methode=calib_methode,
        vizier_dict=vizier_dict,
        file_calib=file_calib,
        mag_range=mag_range,
    )
    pipeline_config.verbose = verbose

    run_transformation_extraction_pipeline(
        observation,
        filter_list,
        img_dirs,
        outdir,
        pipeline_config,
    )

    write_field_transformation_table(
        observation,
        nameobj=nameobj,
        filter_list=filter_list,
        valid_calibs=valid_calibs,
        outdir=outdir,
        extraction_config=pipeline_config,
        weights=weights,
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
