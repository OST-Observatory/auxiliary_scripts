#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
    Bin all images from a specific directory
"""

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
image_path: str = '?'

#   Output directory
output_directory: str = '?'

#   Binning factor
binning: int = 2

############################################################################
####                            Libraries                               ####
############################################################################

import time

from pathlib import Path

import ccdproc as ccdp

from astropy.nddata import block_reduce

from ost_photometry import checks, terminal_output

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    start_time = time.time()
    terminal_output.print_header("Bin images")

    #   Check input and output directory
    file_path = checks.check_pathlib_path(image_path)
    checks.check_output_directories(output_directory)
    out_path = Path(output_directory)

    #   Get image file collection
    ifc = ccdp.ImageFileCollection(file_path)

    #   Loop over files
    bar = terminal_output.progress_bar(len(ifc.files), desc="Bin")
    for ccd, file_name in ifc.ccds(
        ccd_kwargs={'unit': 'adu'},
        return_fname=True,
        ):

        #   Bin image
        #ccd = ccdp.transform_image(
            #ccd,
            #block_reduce,
            #block_size = binning,
            #func=np.mean,
            #)
        ccd = ccdp.block_average(ccd, binning)

        #   Correct Header
        ccd.meta['XBINNING'] = binning
        ccd.meta['YBINNING'] = binning
        ccd.meta['EXPTIME'] = ccd.meta['EXPTIME']/binning/binning
        ccd.meta['EXPOSURE'] = ccd.meta['EXPOSURE']/binning/binning
        ccd.meta['INFO_0'] = 'Software binned using numpy mean function'
        ccd.meta['INFO_1'] = '    Exposure time scaled accordingly'

        #   Save the result
        ccd.write(out_path / file_name, overwrite=True)

        bar.update(1)
    bar.close()

    terminal_output.print_done(start_time, {"Output": out_path})


