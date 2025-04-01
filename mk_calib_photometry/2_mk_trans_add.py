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
#calib_methode = 'URAT1       # URAT1 Catalog (Zacharias+ 2015)
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
#   Critical radius outside which correlations are rejected
dcr     = 3
#dcr     = 5
dcr     = 7

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
#   Create a WCS solution for all images
#mk_all_wcs  = False
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
#   Remove cosmics
#
#   Bool
rmcos = True
rmcos = False

#   Parameters:
objlim  = 5.
sigclip = 4.0

#   Needs to changed in the future
camera = 'QHY600M'


###
#   Camera specific parameters
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


###
#   Expert interface
#   -> add here, if images from more than 2 filters should be reduced
#
filter_list = [filter_1, filter_2]
img_dirs    = {filter_1:img_1, filter_2:img_2}
sigma_psf   = {filter_1:sigma, filter_2:sigma}

############################################################################
####                            Libraries                               ####
############################################################################

import sys

from os.path import join

from shutil import copy

import numpy as np

import time

import warnings
warnings.filterwarnings('ignore')

from photutils.background import MADStdBackgroundRMS

from astropy.stats import SigmaClip
import astropy.units as u
from astropy.table import Table

import multiprocessing as mp

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
    #   Extract flux
    analyze.extract_flux(
        img_container,
        filter_list,
        nameobj,
        img_dirs,
        outdir,
        wcs_method,
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
        rmcos=rmcos,
        objlim=objlim,
        readnoise=readnoise,
        sigclip=sigclip,
        satlevel=satlevel,
        plot_ifi=plot_ifi,
        plot_test=plot_test,
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
            for i in range(0,len(calib_fil)):
                key = calib_fil[i]

                #   Set up filter list
                filt_list    = [calib_fil[0], calib_fil[1]]

                #   Add air mass and object to the calibration table
                tbl_trans['airmass_'+key] = img_container.ensembles[key].get_air_mass()


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


                tbl_trans['jd'] = img_container.ensembles[key].get_obs_time()

    #   Write table and check output directories
    tbl_trans.write(
        outdir+'/tables/trans_para_'+nameobj.replace(' ','_')+'.dat',
        format='ascii',
        overwrite=True,
        )

print(style.bcolors.OKGREEN+"   Done"+style.bcolors.ENDC)
print("--- %s minutes ---" % ((time.time() - start_time)/60.))
