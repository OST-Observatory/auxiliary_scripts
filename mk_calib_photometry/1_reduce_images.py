#! /usr/bin/python
# -*- coding: utf-8 -*-

'''
    First part of the reduction pipeline for data taken within the
    scope of the C7 observation of the astrophysics lab course at
    Potsdam University.

    All files can be given in one directory called 'rawfiles'. Alternatively
    images can be sorted into the following directory structure:
        * Images of the object
        * Dark frames
        * Flatfields
    If they are sorted into directories the FITS header keywords will be
    checked for consistency.

    Images in sub folders will be recognized, but only one level is considered.

   Version
   -------
        0.1   (13.01.2021)
           - adapted from ./1_add_images.py of N2 (see change log there)
        0.11  (03.02.2021)
           - very small style update
        0.12  (18.08.2020)
           - small update ans api adjustments
        0.3  (10.02.2022)
           - complete rewrite using ccdproc
'''

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#######################  Simple folder structure  ##########################
rawfiles = '?'

##########################  Individual folders  ############################
### Path to the bias -- If set to '?', bias exposures are not used.
bias = '?'

### Path to the darks
darks = '?'

### Path to the flats
flats = '?'

### Path to the images
imgs  = '?'


############################################################################
####             Additional options: only edit if necessary             ####
############################################################################

#   Dictionary with file type infos (use None for default, or dict from calibration_parameters.get_image_types())
#   For reduce_main: dict with keys bias/dark/flat/light, values are lists of FITS imagetyp strings
img_type = None  # Uses default from calibration_parameters

#  Path to store the output (will usually be 'output',
#  but it can be changed as needed).
outdir='output/'

##   Verbose output
#verbose = True
verbose = False

#   Number of cores used for multiprocessing (None = half the CPUs)
n_cores_multiprocessing = None


###
#   Remove cosmics
#
#   Bool:
rmcos = True
#rmcos = False

#   Parameters:
objlim  = 5.
sigclip = 5.0
sigclip = 4.0


###
#   Stack images
#
stack = True
stack = False


###
#   Frame quality
#
#   Measure FWHM, roundness, number of stars and sky background of every
#   reduced frame (output/frame_quality.ecsv, FITS headers, QC plots).
measure_frame_quality = True

#   Reject poor frames before the alignment; None keeps every frame.
#   Example: drop frames with a strongly deviating FWHM or almost no stars.
frame_selection = None
#frame_selection = {"fwhm_sigma_clip": 3.0, "n_stars_min": 20}

#   Alignment reference: "best_fwhm" (sharpest frame) or "first".
reference_image_selection = "best_fwhm"


###
#   Camera specific parameters
#
#   Read noise, system gain, dark current and saturation level are taken
#   from the FITS header and the bundled camera catalog
#   (ost_photometry/data/cameras.json). Set a value here only to override.
gain      = None
readnoise = None
dark_rate = None
satlevel  = None


############################################################################
####                            Libraries                               ####
############################################################################

import tempfile

import warnings
warnings.filterwarnings('ignore')

from ost_photometry.reduce import redu, utilities

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    ###
    #   Prepare directories and make checks
    #
    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    #   Prepare directories
    rawfiles = utilities.prepare_reduction(
        outdir,
        bias,
        darks,
        flats,
        imgs,
        rawfiles,
        temp_dir,
        image_type=img_type,
    )

    ###
    #   Reduce images
    #
    redu.reduce_main(
        rawfiles,
        outdir,
        image_type_dir=img_type,
        gain=gain,
        read_noise=readnoise,
        dark_rate=dark_rate,
        rm_cosmic_rays=rmcos,
        saturation_level=satlevel,
        limiting_contrast_rm_cosmic_rays=objlim,
        sigma_clipping_value_rm_cosmic_rays=sigclip,
        debug=verbose,
        stack_images=stack,
        n_cores_multiprocessing=n_cores_multiprocessing,
        measure_frame_quality=measure_frame_quality,
        frame_selection=frame_selection,
        reference_image_selection=reference_image_selection,
    )
