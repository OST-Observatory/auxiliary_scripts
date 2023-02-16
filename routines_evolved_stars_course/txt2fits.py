#!/usr/bin/env python3

import os, glob, argparse
import pandas as pd
import numpy as np
import fitsio

description = '''\
description:
  convert all ascii spectra with some suffix to fuse_register readable .fits
  columns are assumed to be: wave flux [error] [quality]
    '''
epilog = '''\
examples:
  txt2fits.py ./ '*.txt'
  txt2fits.py ./ '*.spek' -err -q
  txt2fits.py ./ '*.spek' -err -q -mq
    '''
parser = argparse.ArgumentParser(
          prog='txt2fits.py',
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description=description,
          epilog=epilog)

parser.add_argument("-err","--errors",
                    action='store_true', default=False,
                    help="Also use error colum")
parser.add_argument("-q","--quality",
                    action='store_true', default=False,
                    help="Also use quality colum")
parser.add_argument("-mq","--mask_quality",
                    action='store_true', default=False,
                    help="rm rows with quality < 0")
parser.add_argument("--multiply",
                    nargs='?', default=False,
                    help="Multiply flux (and error) columns by constant")
parser.add_argument('argv', nargs='*')
args = parser.parse_args()

path_to_structure   =   args.argv[0]
suffix              =   args.argv[1]
dirname  =  os.path.dirname(path_to_structure)
#print(suffix)
#print(dirname)

def get_filenames(dirname,suffix):
    filenames =   []
    searchname = dirname+'/'+suffix
    print(searchname)
    for filename in glob.iglob(searchname):
        #print(filename)
        #if os.path.basename(filename).startswith('HI'):
        filenames.append(filename)
    return(filenames)

def read_to_df(filename):
    df  = pd.read_table(filename,
                     sep = "\s+",
                     header=None,
                     skiprows=0,
                     comment='#',
                     dtype=float,
                     float_precision='high')
    
    # rename to names used in fuse_register
    # WAVE, FLUX, ERROR, COUNTS, WEIGHTS, BKGD, QUALITY
    if args.errors:
        rename_dict = {0: 'WAVE',
                       1: 'FLUX',
                       2: 'ERROR'}
        if args.quality:
            rename_dict = {0: 'WAVE',
                           1: 'FLUX',
                           2: 'ERROR',
                           3: 'QUALITY'}
    else:
        if args.quality:
            rename_dict = {0: 'WAVE',
                           1: 'FLUX',
                           2: 'QUALITY'}
        else:
            rename_dict = {0: 'WAVE',
                           1: 'FLUX'}

    df = df.rename(rename_dict, axis='columns')

    if args.multiply:
        df['FLUX'] = df['FLUX']*float(args.multiply)
        if args.errors:
            df['ERROR'] = df['ERROR']*float(args.multiply)

    if args.quality:
        df.QUALITY = df.QUALITY.astype(np.int64)

    # apply filters
    if args.mask_quality:
        df = df[df['QUALITY'] >= 0]

        # QUALITY is between 0 and 100 in fuse_register (and FUSE data)
        df['QUALITY'] = 100

    #df = df[np.isfinite(df['FLUX']) == True]
    #print(df)

    return df

def df_to_fits(dataframe,name):
    result_file = fitsio.FITS(name, 'rw')
    result_file.write(dataframe.to_records(index=False))
    result_file.close()

def df_to_txt(dataframe,name):
    with open(name, 'w+') as file:
        file.write(dataframe.to_string(index=False))

filenames = get_filenames(dirname,suffix)
for item in filenames:
    print(item)
    fitsname    = item.rsplit('.', 1)[0] + '.fits'
    print(fitsname)
    df          = read_to_df(item)
    print(df[:2])
    df_to_fits(df,fitsname)
