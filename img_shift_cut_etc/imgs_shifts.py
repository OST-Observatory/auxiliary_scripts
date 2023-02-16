#! /usr/bin/python3
# -*- coding: utf-8 -*-

'''
    Load images of specific formats and calculate shifts to a reference
    image
'''


############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

#   Path to the images
path_in  = '?/'

#   Output directory
path_out = '?'

#   Allowed input file formats
#formats = [".tiff", ".TIFF"]
#formats = [".CR2"]
formats = ["?"]

#   Reference image
ref_id = 0


###
#   Image mask -> used to mask image areas that could spoil the
#                 cross correlation such as the moon during a
#                 solar eclipse
#   Mask image?
bool_mask = True

#   Points to define the mask -> the area enclosed by the points
#                                will be masked
mask_points = [[3000, 0],[3000, 5183]]

#   Add a point for the upper right corner
upper_right = False
#   Add a point for the lower right corner
lower_right = True
#   Add a point for the lower left corner
lower_left  = True


###
#   Apply a heaviside function to the image
bool_heavy = False
bool_heavy = True


###
#   Additional cuts to the images
#
#   Upper edge
ys_cut = 0
#   Lower edge
ye_cut = 0
#   Left edge
xs_cut = 0
#   Right edge
xe_cut = 0


###
#   Plot options
#
#   Plot the image mask and reference image
plot_mask = True

#   Plot cut images
plot_cut  = True
#   ID of the image to plot
id_img    = 10


############################################################################
####                            Libraries                               ####
############################################################################

import sys
import os

import numpy as np
import matplotlib.pyplot as plt

from skimage import data
from skimage.registration import phase_cross_correlation
from skimage.draw import polygon2mask
from skimage.io import (imread, imread_collection, imsave, imshow)

from ost import checks, aux_base
from ost.reduce import aux

#from PIL import Image
#from rawkit.raw import Raw
import rawpy

############################################################################
####                               Main                                 ####
############################################################################

#   Check if output directory exists
checks.check_out(path_out)

#   Make file list
sys.stdout.write("\rRead images...\n")
fileList, nfiles = aux_base.mkfilelist(path_in, formats=formats, addpath=True)

#print(fileList)

#   Read images
if ".CR2" in formats:
    im = []
    for i, filename in enumerate(fileList):
        id_c = i+1
        sys.stdout.write("\rImage %i/%i" % (id_c, nfiles))
        sys.stdout.flush()

        buffered_image = rawpy.imread(filename)
        rgb = buffered_image.postprocess()
        im.append(np.array(rgb))
        #raw_image = Raw(filename)
        #buffered_image = np.array(raw_image.to_buffer())
        #image = Image.frombytes(
            #'RGB',
            #(raw_image.metadata.width, raw_image.metadata.height),
            #buffered_image
            #)
        #im.append(numpy.array(image))

    sys.stdout.write("\n")

else:
    im = imread_collection(fileList)

#print(im)
print(im[0].shape)

#   Image shape
image_shape = im[0].shape
nx = image_shape[1]
ny = image_shape[0]

#   Bit depth
if im[0].dtype == 'uint8':
    bit_depth = 8
elif im[0].dtype == 'uint16':
    bit_depth = 16
else:
    print('Caution: Bit depth could not be determined.')
    print('Data type is ', im[0].dtype)
    print('Bit depth set to 8')
    bit_depth = 8


###
#   Calculate image mask -> the masked area will not be considered,
#                           while calculating the shifts
mask = np.ones((ny,nx), dtype=bool)
if bool_mask:
    sys.stdout.write("\rCalculate image mask...\n")
    if upper_right:
        mask_points.append([0,nx])
    if lower_right:
        mask_points.append([ny,nx])
    if lower_left:
        mask_points.append([ny,0])

    #   Calculate mask from the specified points
    mask = np.invert(polygon2mask((ny,nx), mask_points))


###
#   Calculate shifts
#
#   Prepare an array for the shifts
shifts = np.zeros((2,nfiles), dtype=int)

#   Loop over number of files
sys.stdout.write("\rCalculate shifts...\n")
for i in range(0, nfiles):
    #   Write current status
    id_c = i+1
    sys.stdout.write("\rImage %i/%i" % (id_c, nfiles))
    sys.stdout.flush()

    #   "Normalize" image & calculate heaviside function for the image
    if bool_heavy:
        #reff = np.heaviside(im[ref_id][:,:,0]/255, 0.03)
        #test = np.heaviside(im[i][:,:,0]/255, 0.03)
        #reff = np.heaviside(im[ref_id][:,:,0]/2**(bit_depth), 0.03)
        #test = np.heaviside(im[i][:,:,0]/2**(bit_depth), 0.03)
        reff = np.heaviside(im[ref_id][:,:,0]/2**(bit_depth), 0.01)
        test = np.heaviside(im[i][:,:,0]/2**(bit_depth), 0.01)
    else:
        reff = im[ref_id][:,:,0]
        test = im[i][:,:,0]

    #   Calculate shifts
    shifts[:,i] = phase_cross_correlation(reff, test, reference_mask=mask)

    #   Plot reference image and image mask
    if i == 0 and plot_mask:
        fig = plt.figure(figsize=(12, 7))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1)

        ax1.imshow(reff, cmap='gray')
        ax1.set_axis_off()
        ax1.set_title('Reference image')

        ax2.imshow(mask, cmap='gray')
        ax2.set_axis_off()
        ax2.set_title('Image mask')

        plt.show()


sys.stdout.write("\n")

#   Ensure shifts are Integer and that they have the right sign
shifts = shifts.astype(int)*-1


###
#   Cut pictures and add them
#
img_cut = {}
for i in range(0,nfiles):
    #   Write status to console
    id_c = i+1
    sys.stdout.write("\rApply shift to image %i/%i" % (id_c, nfiles))
    sys.stdout.flush()

    xs, xe, ys, ye =  aux.make_index_from_shifts(
        shifts,
        i,
        )

    #   Actual image cutting
    img_cut[i] = im[i][ys-ys_cut:ye-ye_cut, xs-xs_cut:xe-xe_cut]

    #   Plot reference and offset image
    if i == id_img and plot_cut and ref_id <= id_img:
        fig = plt.figure(figsize=(12, 7))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1)

        ax1.imshow(img_cut[ref_id])
        ax1.set_axis_off()
        ax1.set_title('Reference image')

        ax2.imshow(img_cut[id_img])
        ax2.set_axis_off()
        ax2.set_title('Offset image')

        plt.show()

    #   Write image
    new_name = 'shift_'+os.path.basename(fileList[i])
    imsave(os.path.join(path_out,new_name), img_cut[i])

    #new_name = 'shift_'+basename(fileList[i]).split('.')[0]
    #imsave(join('output/jpeg',new_name)+'.jpg', img_cut[i])


sys.stdout.write("\n")







