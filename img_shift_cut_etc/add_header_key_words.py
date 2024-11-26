#! /home/pollux/.virtualenvs/photo/bin/python
# -*- coding: utf-8 -*-

"""
    Add header keywords to image
"""

from astropy.nddata import CCDData

#------------------------------------------------------------------------------

path: str = '?'

#------------------------------------------------------------------------------

def get_basename(path):
    """
        Determine base name without ending from a file path. Accounts for
        multiple dots in the file name.

        Parameters
        ----------
        path            : `string` or `pathlib.Path` object
            Path to the file

        Returns
        -------
        base_name        : `string`
            Basename without ending
    """
    name_parts = str(path).split('/')[-1].split('.')[0:-1]
    if len(name_parts) == 1:
        base_name = name_parts[0]
    else:
        base_name = name_parts[0]
        for part in name_parts[1:]:
            base_name = base_name+'.'+part

    return base_name

#------------------------------------------------------------------------------
#   Read image
image  = CCDData.read(path)

#   Get observation dat
date   = image.meta['DATE-OBS']

#   Get target
target = image.meta['OBJECT']
target = 'GRB220412A'

#   Get gain
gain   = image.meta['EGAIN']

#   Add keywords to Header
image.meta['USERNAME']      = 'OST'
image.meta['INSTRU']        = 'CDK'
image.meta['OBSDATE']       = date
#image.meta['FILTER system'] = 'bessell'
image.meta['FILTER-S']      = 'bessell'
image.meta['TARGET']        = target
image.meta['STACK']         = 1
image.meta['GAIN']          = gain

#   Set exposure time
#image.meta['EXPTIME']       = 56*180
#image.meta['EXPOSURE']      = image.meta['EXPTIME']

#   Set file name
basename  = get_basename(path)
file_name = basename+'_mod.fit'

#   Write file to disk
image.write(file_name, overwrite=True)


