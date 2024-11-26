#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
    Trim all images from a specific directory
"""

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
image_path: str = '?'

#   Output directory
output_directory: str = '?'

############################################################################
####                            Libraries                               ####
############################################################################

import sys

from pathlib import Path

import ccdproc as ccdp

from ost_photometry import checks

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Check input and output directory
    file_path = checks.check_pathlib_path(image_path)
    checks.check_output_directories(output_directory)
    out_path = Path(output_directory)

    #   Get image file collection
    ifc = ccdp.ImageFileCollection(file_path)

    #   Loop over files
    i = 0
    for img, file_name in ifc.ccds(
        ccd_kwargs={'unit': 'adu'},
        return_fname=True,
        ):

        #   Trim image
        img = ccdp.trim_image(img[300:1840, 400:-400])

        #   Save the result
        img.write(out_path / file_name, overwrite=True)

        #   Write status to console
        i += 1
        sys.stdout.write("Trim image %i" % i)
        sys.stdout.flush()
    sys.stdout.write("\n")


