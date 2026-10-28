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

import time

from pathlib import Path

import ccdproc as ccdp

from ost_photometry import checks, terminal_output

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    start_time = time.time()
    terminal_output.print_header("Trim images")

    #   Check input and output directory
    file_path = checks.check_pathlib_path(image_path)
    checks.check_output_directories(output_directory)
    out_path = Path(output_directory)

    #   Get image file collection
    ifc = ccdp.ImageFileCollection(file_path)

    #   Loop over files
    with terminal_output.progress_bar(len(ifc.files), desc="Trim") as bar:
        for img, file_name in ifc.ccds(
            ccd_kwargs={'unit': 'adu'},
            return_fname=True,
            ):

            #   Trim image
            img = ccdp.trim_image(img[300:1840, 400:-400])

            #   Save the result
            img.write(out_path / file_name, overwrite=True)
            bar.update(1)

    terminal_output.print_done(start_time, {"Output": out_path})


