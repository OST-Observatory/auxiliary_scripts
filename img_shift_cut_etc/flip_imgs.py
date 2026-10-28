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

import os

import time

import numpy as np

from tifffile import (imread, imsave)

from ost_photometry import checks, terminal_output, utilities

############################################################################
####                               Main                                 ####
############################################################################

start_time = time.time()
terminal_output.print_header("Flip images")

#   Check if output directory exists
checks.check_output_directories(output_directory)

#   Make file list
terminal_output.print_to_terminal("Read images...")
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
bar = terminal_output.progress_bar(n_files, desc="Flip")
for i in range(0, n_files):
    #   Extract file name
    new_name = os.path.basename(fileList[i])

    #   Write images
    imsave(os.path.join('flip_imgs_test',new_name), im[i])

    bar.update(1)
bar.close()

terminal_output.print_done(start_time, {"Output": 'flip_imgs_test'})


