#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
    Bin all images from a specific directory
"""

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
path = '?'

#   Output directory
outdir = '?'

#   Binning factor
binning = 2

############################################################################
####                            Libraries                               ####
############################################################################

import sys

from pathlib import Path

import numpy as np

import ccdproc as ccdp

from astropy.nddata import block_reduce

from ost_photometry import checks

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Check input and output directory
    file_path = checks.check_pathlib_path(path)
    checks.check_output_directories(outdir)
    out_path = Path(outdir)

    #   Get image file collection
    ifc = ccdp.ImageFileCollection(file_path)

    #   Loop over files
    i = 0
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

        #   Write status to console
        i += 1
        sys.stdout.write("\Bin image %i" % i)
        sys.stdout.flush()
    sys.stdout.write("\n")


