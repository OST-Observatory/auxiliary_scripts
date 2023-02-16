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

#   Preserve reference objects in cross correlation between different images
#   -> images without reference objects will be rejected
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

#   Plot sigma clipped magnitudes
#plot_sigma = True
plot_sigma = False


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

from os.path import join

import time

import warnings
warnings.filterwarnings('ignore')

from astropy.table import Table
from ost_photometry import checks
from ost_photometry import style
from ost_photometry.analyze import (
    analyze,
    plot,
    calib,
    trans,
    aux,
    correlate,
    )


############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Set start time
    start_time = time.time()


    ###
    #   Initialize image ensemble container
    #
    img_container = analyze.image_container()


    ###
    #   Check output directories
    #
    checks.check_out(
        outdir,
        join(outdir,'tables'),
        )


    ###
    #   Check image directories
    #
    checks.check_dir(img_dirs)


    #   Outer loop over all filter
    for filt in filter_list:
        print(
            style.bcolors.HEADER
            +"   Analyzing "+filt+" images"
            +style.bcolors.ENDC
            )


        #   Initialize image ensemble object
        img_container.ensembles[filt] = analyze.image_ensemble(
                filt,
                nameobj,
                img_dirs[filt],
                outdir,
                ref_ID,
                )


        ###
        #   Find the WCS solution for the image
        #
        aux.find_wcs(
            img_container.ensembles[filt],
            ref_ID,
            method=wcs_method,
            indent='         ',
            )


        ###
        #   Main extraction of object positions and object fluxes
        #   using multiprocessing
        #
        analyze.extract_multiprocessing(
            img_container.ensembles[filt],
            ncores,
            sigma_bkg,
            sigma_psf,
            multi_start=multi_start,
            size_epsf=size_epsf,
            frac_epsf_stars=frac_epsf_stars,
            oversampling=oversampling,
            maxiters=maxiters,
            methode=methode,
            multi=multi,
            multi_grouper=multi_grouper,
            strict_cleaning=strict_cleaning,
            min_eps_stars=min_eps_stars,
            strict_eps=strict_eps,
            photometry=photometry,
            rstars=rstars,
            rbg_in=rbg_in,
            rbg_out=rbg_out,
            r_unit=r_unit,
            plot_ifi=plot_ifi,
            plot_test=plot_test,
            )


        ###
        #   Correlate results from all images, while preserving the
        #   calibration stars
        #
        analyze.correlate_preserve_calibs(
            img_container.ensembles[filt],
            filter_list,
            calib_methode=calib_methode,
            mag_range=mag_range,
            vizier_dict=vizier_dict,
            calib_file=file_calib,
            dcr=dcr,
            option=option,
            verbose=verbose,
            maxid=maxid,
            ref_ID=ref_ID,
            nmissed=nmissed,
            bfrac=bfrac,
            s_refOBJ=s_refOBJ,
            plot_test=plot_test,
            )


        ###
        #   Make new table to add the median of the flux
        #
        aux.add_median_table(img_container.ensembles[filt])


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
            for i in range(0,len(calib_fil)):
                key = calib_fil[i]

                #   Set up filter list
                filt_list    = [calib_fil[0], calib_fil[1]]

                #   Add air mass and object to the calibration table
                tbl_trans['airmass_'+key] = [
                    img_container.ensembles[key].median_air_mass()
                    ]

                ###
                #   Correlate the results from the different filter and
                #   determine transformation coefficients
                #
                trans.calculate_trans(
                    img_container,
                    key,
                    filt_list,
                    tbl_trans,
                    weights=weights,
                    dcr=dcr,
                    option=option,
                    calib_methode=calib_methode,
                    vizier_dict=vizier_dict,
                    calib_file=file_calib,
                    mag_range=mag_range,
                    )

                tbl_trans['jd'] = [
                    img_container.ensembles[key].median_obs_time()
                    ]


    #   Write table and check output directories
    tbl_trans.write(
        outdir+'/tables/trans_para_'+nameobj.replace(' ','_')+'.dat',
       format='ascii',
       overwrite=True,
       )

print(style.bcolors.OKGREEN+"   Done"+style.bcolors.ENDC)
print("--- %s minutes ---" % ((time.time() - start_time)/60.))
