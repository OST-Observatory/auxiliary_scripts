#! /usr/bin/python3
# -*- coding: utf-8 -*-

'''
    Load images of specific formats and calculate shifts to a reference
    image
'''


############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
#   The path can be a given as a string to a single directory or a list of
#   strings to multiple directories.
path_in = [
    'out_indi_record15-07-331/video_imgs/',
    'out_indi_record15-16-242/video_imgs/',
    'out_indi_record15-35-093/video_imgs/',
    'out_indi_record15-40-334/video_imgs/',
    ]

#   Output directory
path_out = 'out_all/'

#   Allowed input file formats
formats = [".tiff", ".TIFF"]
# formats = [".fit", ".fits"]

#   Output format
out_format = ".tiff"
#out_format = ".fit"

#   ID of the reference image
ref_id = 0


###
#   Extend or cut image -> The images can be either cut to the common
#                          field of view or extended to a field of view
#                          in which all images fit
#                       -> Limitation: ``extend`` is currently not available
#                                      for FITS files
#
mode = 'extend'
mode = 'cut'


###
#   Image mask -> used to mask image areas that could spoil the
#                 cross correlation such as the moon during a
#                 solar eclipse
#
#   Mask image?
bool_mask = True
bool_mask = False

#   Points to define the mask -> the area enclosed by the points
#                                will be masked
mask_points = [[0, 0]]

#   Add a point for the upper right corner
upper_right = True
#   Add a point for the lower right corner
lower_right = True
#   Add a point for the lower left corner
lower_left  = False


###
#   Apply a heavyside function to the image
#
bool_heavy = False
bool_heavy = True

#   Background offset
offset = 2500
offset = 6800


###
#   Additional cuts to the images
#
#   Upper edge
ys_cut = 0
#   Lower edge
ye_cut = 0
#   Left edge
xs_cut = 0
#   Right edge
xe_cut = 0


###
#   Plot options
#
#   Plot the image mask and reference image
plot_mask = True
#plot_mask = False

#   Plot cut images
plot_cut  = True
plot_cut  = False
#   ID of the image to plot
id_img    = 10


############################################################################
####                            Libraries                               ####
############################################################################

import sys
import os

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from skimage import data
from skimage.registration import phase_cross_correlation
from skimage.draw import polygon2mask
from skimage.io import imread, imread_collection, imsave, imshow

import checks
import aux

import registration


############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Check if command line arguments are given
    #   -> Used to overwrite path and format variables from above
    if len(sys.argv) == 5:
        path_in = sys.argv[1]
        path_out = sys.argv[2]
        formats = sys.argv[3]
        out_format = sys.argv[4]

    #   Check if the path/s to the images are given as a string or list
    #   In case of a string replace it with a list
    if isinstance(path_in, str):
        path_in = [path_in]
    elif not isinstance(path_in, list):
        raise TypeError('The path to the images ({}) is neither a list nor a'
                        'string => Exit'.format(path_in))
    #   The same check for the input formats
    if isinstance(formats, str):
        formats = [formats]
    elif not isinstance(formats, list):
        raise TypeError('The input formats ({}) are specified neither as a list '
                        'nor as a string => Exit'.format(formats))

    #   Calculate and apply shifts
    if '.fit' in formats or '.fits' in formats:
        registration.cal_img_shifts_fits(
            path_in,
            path_out,
            formats,
            out_format,
            ref_id=ref_id,
            mode=mode,
            bool_mask=bool_mask,
            mask_points=mask_points,
            upper_right=upper_right,
            lower_right=lower_right,
            lower_left=lower_left,
            bool_heavy=bool_heavy,
            offset=offset,
            ys_cut=ys_cut,
            ye_cut=ye_cut,
            xs_cut=xs_cut,
            xe_cut=xe_cut,
            plot_mask=plot_mask,
            plot_cut=plot_cut,
            id_img=id_img,
            )
    else:
        registration.cal_img_shifts_normal(
            path_in,
            path_out,
            formats,
            out_format,
            ref_id=ref_id,
            mode=mode,
            bool_mask=bool_mask,
            mask_points=mask_points,
            upper_right=upper_right,
            lower_right=lower_right,
            lower_left=lower_left,
            bool_heavy=bool_heavy,
            offset=offset,
            ys_cut=ys_cut,
            ye_cut=ye_cut,
            xs_cut=xs_cut,
            xe_cut=xe_cut,
            plot_mask=plot_mask,
            plot_cut=plot_cut,
            id_img=id_img,
            )
