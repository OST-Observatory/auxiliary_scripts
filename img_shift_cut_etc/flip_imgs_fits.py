#! /usr/bin/python3
# -*- coding: utf-8 -*-

'''
    Flip all images from a specific directory
'''

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
file_path = '?'

#   Output directory
out_path = '?'

############################################################################
####                            Libraries                               ####
############################################################################

import ccdproc as ccdp

from pathlib import Path

from ost.reduce import aux

############################################################################
####                               Main                                 ####
############################################################################

#   Setup image file collections
ifc = ccdp.ImageFileCollection(file_path)

#   Flip images
aux.flip_img(ifc, Path(out_path))
