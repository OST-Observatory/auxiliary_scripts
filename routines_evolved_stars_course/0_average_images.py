#!/usr/bin/env python3

def average_images(flatdark_dir,masterflat_output_name):
    import os
    import numpy as np
    from astropy.io import fits
    flatdark_list  = os.listdir(flatdark_dir)
    with fits.open(flatdark_dir+"/"+flatdark_list[0]) as setup_flatdark_image:
        setup_flatdark_image_data = setup_flatdark_image[0].data
        size_flatdark_x = len(setup_flatdark_image_data)
        size_flatdark_y = len(setup_flatdark_image_data[1])
    flatdark_iarray = np.zeros((size_flatdark_x,size_flatdark_y,len(flatdark_list)))
    i = 0
    for element in flatdark_list:
        with fits.open(flatdark_dir+"/"+element) as flatdark_image:
            flatdark_image_data = flatdark_image[0].data#.astype(int)
            flatHeader = flatdark_image[0].header
            flatdark_iarray[...,flatdark_list.index(element)] = flatdark_image_data
            i = i+1
            print("Reading",element)
    flatdark_iarray_av = np.median(flatdark_iarray.astype(int), axis=2)
    #print(flatdark_iarray_av.shape)
    hdu = fits.PrimaryHDU(flatdark_iarray_av,flatHeader)
    hdu.scale('int16', '', bzero=0)
    hduList = fits.HDUList([hdu])
    hduList.writeto(masterflat_output_name, overwrite=True)
    print("Saved to",masterflat_output_name)

def average_images_2(flatdark_dir,masterflat_output_name):
    from netCDF4 import Dataset
    import numpy as np
    from PIL import Image
    from astropy.io import fits
    import glob
    import os

    images = glob.glob(flatdark_dir+'*.fits')

    with fits.open(images[0]) as setup_flatdark_image:
        setup_flatdark_image_data = setup_flatdark_image[0].data
        size_flatdark_x = len(setup_flatdark_image_data)
        size_flatdark_y = len(setup_flatdark_image_data[1])

    WIDTH = size_flatdark_x
    HEIGHT = size_flatdark_y
    root = Dataset('test.nc', 'w')
    root.createDimension('x', WIDTH)
    root.createDimension('y', HEIGHT)
    root.createDimension('frames', None)

    # 'u2' -> 16-bit unsigned integer
    img = root.createVariable('image', 'u2', ('frames','x','y'),fill_value=65535)

    for i,element in enumerate(images):
        with fits.open(element) as flatdark_image:
            flatdark_image_data = flatdark_image[0].data#.astype(int)
            flatHeader = flatdark_image[0].header
            im_array  = np.asarray(flatdark_image_data)
            img[i] = im_array
            i=i+1
            print(i)

    median = np.median(img, axis=0)
    os.system('rm test.nc')

    im = Image.fromarray(np.uint16(median))
    hdu = fits.PrimaryHDU(median,flatHeader)
    hdu.scale('int16', '', bzero=0)
    hduList = fits.HDUList([hdu])
    hduList.writeto(masterflat_output_name, overwrite=True)

def main():
    import argparse

    description = '''\
    description:
      average all .fits images in the specified directory (mean is the default)
      all images must have the same dimensions
        '''
    epilog = '''\
    examples:
      0_average_images.py path_to_dir/
      0_average_images.py hgar/ -o hgar_stacked.fit
      0_average_images.py science_1200s/ -o science_1200s_stacked.fit
        '''

    parser = argparse.ArgumentParser(
          prog='average_images.py',
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description=description,
          epilog=epilog)
    parser.add_argument('-o',"--outname",
                        nargs=1, default=['stacked.fit'],
                        help="Specify savename")
    parser.add_argument('-m',"--mode",
                        nargs=1, default=['median'],
                        help="image reduction mode: ['mean','median']")
    parser.add_argument('argv', nargs='*')
    args = parser.parse_args()

    dir_stack = args.argv[0]
    outname = args.outname[0]
    average_images(dir_stack,outname)
    #average_images_2(dir_stack,outname)

if __name__ == "__main__":
    main()
