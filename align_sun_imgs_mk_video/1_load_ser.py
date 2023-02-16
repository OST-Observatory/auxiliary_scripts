#! /usr/bin/python3
# -*- coding: utf-8 -*-

'''
    Process .ser files
'''

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to video
file_path = '?.ser'

#   Output directory
path_out = 'out/'

#   Image output format
out_format = '.tiff'
# out_format = '.fit'

############################################################################
####                            Libraries                               ####
############################################################################

import os
import sys

from pathlib import Path, PurePath

import planetary_system_stacker as pss
sys.path.append(pss.__path__[0])
from configuration import Configuration
from frames import Frames
import ser_parser

import cv2

import checks
import aux

from skimage.io import imsave

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Check directories
    sys.stdout.write("\rCheck directories...\n")
    path_in  = checks.check_pathlib_Path(file_path)
    checks.check_out(path_out)

    #   Path to the postprocess images
    imgs_path = Path(path_out)
    imgs_path = imgs_path / 'video_imgs'
    imgs_path.mkdir(exist_ok=True)

    #   Load configuration for the Planetary System Stacker
    configuration = Configuration()
    configuration.initialize_configuration(read_from_file=False)

    #   Create the VideoCapture object.
    cap = ser_parser.SERParser(file_path, SER_16bit_shift_correction=True)

    #   Get number of frames
    frame_count = cap.frame_count

    #	Get digits of the number of frames -> for the file name
    digits = len(str(frame_count))

    #	Change the Debayer pattern (is somehow needed in our case)
    cap.header['DebayerPattern'] = cv2.COLOR_BayerBG2BGR

    #	Read all frames
    data = cap.read_all_frames()

    #   Loop over all frames
    for j, img in enumerate(data):
        #   Create new name and path
        new_name = aux.get_basename(file_path)+'_'+str(j+1).rjust(digits, '0')
        out_path = imgs_path / str(new_name+out_format)
        #   Write frame as image
        imsave(
            out_path,
            img,
            )
