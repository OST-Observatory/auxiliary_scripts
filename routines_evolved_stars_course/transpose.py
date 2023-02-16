#!/usr/bin/env python3

import numpy as np
import argparse

def save_arrays_as_ascii(savename,array_list):
    array_list = [np.array(i).astype(float) for i in array_list]
    spec_save = np.transpose([array_list[0],array_list[1]])
    with open(savename,"w") as savefile:
        np.savetxt(savename,spec_save)
        print("Saved to",savename)

description = '''\
description:
  convert ascii file with two rows and many columns
  to ascii file with two columns and many rows
    '''
epilog = '''\
examples:
  transpose.py spectrum_science_sky_f.dat
    '''

parser = argparse.ArgumentParser(
          prog='transpose.py',
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description=description,
          epilog=epilog)
parser.add_argument('argv', nargs='*')
args = parser.parse_args()

readname = args.argv[0]

a, b = np.loadtxt(readname)

savename = readname + '.T'
save_arrays_as_ascii(savename,[a,b])
