#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
    Flip all images from a specific directory
"""

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
file_path: str = '?'

#   Output directory
out_path: str = '?'

############################################################################
####                            Libraries                               ####
############################################################################

import ccdproc as ccdp

from pathlib import Path

from ost_photometry.reduce import utilities

############################################################################
####                               Main                                 ####
############################################################################

#   Setup image file collections
ifc = ccdp.ImageFileCollection(file_path)

#   Flip images
utilities.flip_image(ifc, Path(out_path))
