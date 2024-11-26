#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
    Flip all images from a specific directory
"""

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
image_path: str = '?'

#   Output directory
output_directory: str = '?'

#   Allowed input file formats
formats: list[str] = [".tiff", ".TIFF"]

############################################################################
####                            Libraries                               ####
############################################################################

import sys

import os

import numpy as np

from tifffile import (imread, imsave)

from ost_photometry import checks, utilities

############################################################################
####                               Main                                 ####
############################################################################

#   Check if output directory exists
checks.check_output_directories(output_directory)

#   Make file list
sys.stdout.write("\rRead images...\n")
fileList, n_files = utilities.mk_file_list(
    image_path,
    formats=formats,
    add_path_to_file_names=True,
)

#   Read images
im = imread(fileList, key=0)

#   Flip images
im = np.flip(im, axis=1)
im = np.flip(im, axis=2)

#   Write images to the output directory
for i in range(0, n_files):
    #   Extract file name
    new_name = os.path.basename(fileList[i])

    #   Write images
    imsave(os.path.join('flip_imgs_test',new_name), im[i])

    #   Write status to console
    id_img = i+1
    sys.stdout.write("\rFlip image %i" % id_img)
    sys.stdout.flush()
sys.stdout.write("\n")


