#!/usr/bin/env python3

from astropy.io import fits
import argparse, os, glob
import numpy as np

def get_filenames(dirname,suffix):
    filenames =   []
    searchname = dirname+'/'+suffix
    print(searchname)
    for filename in glob.iglob(searchname):
        filenames.append(filename)
    return(filenames)

def rename_file(fname,key):
    with fits.open(fname) as fimage:
        fheader = fimage[0].header
        exptime = str(fheader[key])

    newname = fname.rsplit(".",1)[0] + "_" + exptime  + "." + fname.rsplit(".",1)[1]

    os.system("cp " + fname + " " + newname)

def rotate_frame(fname):
    with fits.open(fname) as fimage:
        iheader = fimage[0].header
        idata = fimage[0].data
        idata_90 = np.rot90(idata)

    newname = fname.rsplit(".",1)[0] + "_" + "r90"  + "." + fname.rsplit(".",1)[1]
    newname = fname
    hdu = fits.PrimaryHDU(idata_90,iheader)
    hdu.scale('int16', '', bzero=0)
    hduList = fits.HDUList([hdu])
    hduList.writeto(newname, output_verify='exception', overwrite=True)

def main():

    description = "renames fits files to include EXPTIME header"
    epilog = "~/spectroscopic_analysis/scripts/rename_fits.py . '*.fits'"
    parser = argparse.ArgumentParser(
          prog='rename_fits.py',
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description=description,
          epilog=epilog)
    parser.add_argument('argv', nargs='*')
    args = parser.parse_args()

    dirname = args.argv[0]
    suffix = args.argv[1]

    filenames = get_filenames(dirname,suffix)

    for item in filenames:
        #rename_file(item,'EXPTIME')
        rotate_frame(item)

main()
