#! /usr/bin/python3
# -*- coding: utf-8 -*-

'''
    Make video from individual images
'''

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

###
#   Paths and formats definitions
#

#   Path to the images
path_in = 'out_all/cut/'

#   Output directory
path_out = 'out_all/'

#   Allowed input file formats
#formats = [".jpg", ".jpeg", ".JPG", ".JPEG"]
#formats = [".FIT", ".fit", ".FITS", ".fits"]
formats = [".tif", ".tiff", ".TIF", ".TIFF"]

#   Image output format
#out_format = '.jpg'
out_format = '.tiff'


###
#   Trim images
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
#   Find the best images by means of Planetary System Stacker
#
best_img = False
#best_img = True

#   Number of best images to br returned
nimgs  = 1

#   The ''nimgs'' need to be in an interval of size ''window''
window = 1

#   Step size to be evaluated, e.g., ''nimgs=1'' and ''step=20'' means:
#   select the best image out of every 20 images
step   = 5
#step   = 10
step   = 20

#   Add last step even if it is < "step"
add_last_best = False


###
#   Stack best images
#
stack = False
stack = True

#   % of images to be stacked
stack_percent = 20

#   Interval to stack
#stack_interval = 100
#stack_interval = 40
stack_interval = 20
#stack_interval = 10
#stack_interval = 5

#   Drizzling: Interpolate input frames by a drizzle factor
#              Possible values: 1.5x, 2x, 3x, Off
#              A "drizzle" times larger image will be returned.
#              To account for this the images trim values from
#              above will be scaled with "drizzle"
drizzle = 'Off'
#drizzle = '1.5x'

#   The following parameters usually don't need adjustment
#   Noise level (add Gaussian blur) - Range: 0-11
noise = 5
#   Alignment point box width (pixels) for multipoint alignment
box_width = 52
#   Alignment point search width (pixels)
search_width = 20
#   Minimum structure for multipoint alignment
min_struct = 0.07
#   Minimum brightness for multipoint alignment
min_bright = 50
#   Add last step even if it is < "stack_interval"
add_last_stack = False



###
#   Make RGB from Gray scale? (only for not FITS file formats)
#
mkRGB = True

#   Scaling parameters
r_scale = 0.9
g_scale = 0.5
b_scale = 0.5


###
#   Image postprocessing
#
#   Global histogram equalization
global_histo_equalization = True
global_histo_equalization = False

#   Local histogram equalizapath_in  = 'tion
# local_histo_equalization = True
local_histo_equalization = False
#   Footprint for local histogram equalization (disk morphology is used)
disk_size = 30

#   Adaptive histogram equalization
adaptive_histo_equalization = True
adaptive_histo_equalization = False
#   Clip limit for adaptive histogram equalization
clip_limit = 0.03
clip_limit = 0.002

#   log contrast adjustment
log_cont_adjust = False
#   gain for log contrast adjustment
log_gain = 1.

#   Gamma contrast adjustment
gamma_adjust = True
# gamma_adjust = False
#   Gamma value
# gamma = 0.3
# gamma = 0.4
gamma = 0.5
# gamma = 0.6
# gamma = 0.8
# gamma = 1.3
# gamma = 3.

#   Contrast stretching
contrast_stretching = True
contrast_stretching = False
#   Upper percentile for contrast stretching
# upper_percentile = 90.
# upper_percentile = 98.
# upper_percentile = 99.
upper_percentile = 100.
#   Lower percentile for contrast stretching
# lower_percentile = 10.
# lower_percentile = 5.
lower_percentile = 2.
# lower_percentile = 0.

#   Image normalization with ImageMagick
norm_image_magick = True
# norm_image_magick = False

#   Define parameters for postprocessing/sharpening
#   (multiple "layers" are possible)
#
#   Parameters
#   ----------
#       radius          : `float`
#           Radius (in pixels) of the Gaussian sharpening kernel.
#       amount          : `float`
#           Amount of sharpening for this layer.
#       bi_fraction     : `float`
#           Fraction of bilateral vs. Gaussian filter (0.: only Gaussian,
#           1.: only bilateral).
#       bi_range        : `float`
#           Luminosity range parameter of bilateral filter
#           (0 <= bi_range <= 255).
#       denoise         : `float`
#           Fraction of Gaussian blur to be applied to this layer
#           (0.: No Gaussian blur, 1.: Full filter application).
#       luminance_only  : `boolean`
#           True, if sharpening is to be applied to the luminance
#           channel only.
#
#   Usage: [radius, amount, bi_fraction, bi_range, denoise, luminance_only]

#   Example:
postprocess_layers = [
    [1.9, 6., 0.5, 20, 0.8, False,],
    [2.4, 6., 0.35, 20, 0.72, False],
    [3.9, 2.0, 0.35, 20, 0.72, False],
    ]
#postprocess_layers = [
    #[1.9, 4., 0.5, 20, 0.8, False,],
    #[2.4, 4., 0.35, 20, 0.72, False],
    #[3.9, 1.2, 0.35, 20, 0.72, False],
#    ]
postprocess_layers = [
    [1.9, 12., 0.5, 20, 0.8, False,],
    [2.4, 12., 0.35, 20, 0.72, False],
    [3.9, 4., 0.35, 20, 0.72, False],
    ]
postprocess_layers = [
    [1.9, 18., 0.5, 20, 0.8, False,],
    [2.4, 18., 0.35, 20, 0.72, False],
    [3.9, 6., 0.35, 20, 0.72, False],
    ]
postprocess_layers = [
    [1.9, 25., 0.5, 20, 0.8, False,],
    [2.4, 25., 0.35, 20, 0.72, False],
    [3.9, 8., 0.35, 20, 0.72, False],
    ]

############################################################################
####                            Libraries                               ####
############################################################################

import sys
import os

import subprocess

import tempfile

from pathlib import Path

import numpy as np

import ccdproc as ccdp

from astropy.nddata import CCDData

import warnings
warnings.filterwarnings('ignore')

from astropy import log
log.setLevel('ERROR')

from skimage import data
from skimage.io import imsave, imread_collection

#from skimage import exposure
#from skimage.morphology import disk
#from skimage.filters import rank

import planetary_system_stacker as pss
sys.path.append(pss.__path__[0])
from configuration import Configuration
from frames import Frames
from rank_frames import RankFrames
#from planetary_system_stacker.planetary_system_stacker import (
    #Configuration,
    #Frames,
    #RankFrames,
    #)

import postprocess
import checks
import aux

############################################################################
####                               Main                                 ####
############################################################################

#   Sanity check: Stacking and best image selection should not be
#                 selected simultaniously
if best_img and stack:
    print("Stacking and best image selection should not be",
           "selected simultaenously. Choose one!")
    sys.exit()

#   Check directories
sys.stdout.write("\rCheck directories...\n")
path_in  = checks.check_pathlib_Path(path_in)
checks.check_out(path_out)

#   Path to the postprocess images
post_path = Path(Path(path_out) / 'postpro')
post_path.mkdir(exist_ok=True)

#   FITS?
fitsformats = set([".FIT", ".fit", ".FITS", ".fits"])
bfits = len(set(formats).intersection(fitsformats)) >= 1

#   Get file collection
sys.stdout.write("\rRead images...\n")
if bfits:
    ifc = ccdp.ImageFileCollection(path_in)

    #   Apply filter to the image collection
    #   -> This is necessary so that the path to the image directory is added
    #      to the file names. This is required for `calculate_image_shifts`.
    ifc = ifc.filter(SIMPLE=True)

    #   File list
    files = ifc.files

    #   Number of files
    nfiles = len(files)
else:
    files, nfiles = aux.mkfilelist(
        path_in,
        formats=formats,
        addpath=True,
        sort=True,
        )
    if len(files) == 0:
        raise RuntimeError('No files found -> EXIT')
    ifc = imread_collection(files)

#	Get digits of number of files
digits = len(str(nfiles))

#   Get current directory
pwd = os.path.dirname(os.path.abspath(__file__))

#   Get reference image
if bfits:
    ref_img_name = ifc.files[0]
    ref_img_data = CCDData.read(ref_img_name, unit='adu').data
else:
    ref_img_data = ifc[0]

#   Bit depth
if ref_img_data.dtype == 'uint8':
    bit_depth = 8
elif ref_img_data.dtype == 'uint16':
    bit_depth = 16
else:
    print('Caution: Bit depth could not be determined.')
    print('Data type is ', ref_img_data.dtype)
    print('Bit depth set to 16')
    bit_depth = 16


###
#   Stack best images
#
if stack:
    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()
    temp_dir_stack = tempfile.TemporaryDirectory()

    #   Prepare a list for stacked images
    #path_stacked_imgs = []

    j = 1
    for i, img_path in enumerate(files):

        #   Get image file name
        filename = str(img_path).split('/')[-1]

        #   File type
        filetype = filename.split('.')[-1]

        ##   Get image base name
        #basename = aux.get_basename(img_path)

        ##   Compose file name
        #filename = basename+'.fit'

        #   Create a link to the image in the temporary directory
        os.symlink(pwd+'/'+img_path, temp_dir.name+'/'+filename)

        #   If set is reach, find best images
        if ((i != 0 and i%stack_interval == 0) or
            (i+1 == len(files) and i%stack_interval != 0) and add_last_stack):

            #   Construct command for the Planetary System Stack/tmp/tmp4gxsqr6p_pss.tiffer
            #command = "python3 planetary_system_stacker.py " \
            #+ " --out_format tiff -b 4 -s "+str(int(stack_percent)) \
                    #+ " -a 52 -w 20 --normalize_bright --drizzle 1.5x"
                    #+ " --out_format fits -b 4 -s "+str(int(stack_percent)) \
            command = "PlanetarySystemStacker " \
                    + temp_dir.name \
                    + " --out_format tiff -b 4 -s "+str(int(stack_percent)) \
                    + " -a "+str(box_width)+" -w "+str(search_width) \
                    + " --normalize_bright" \
                    + " --align_min_struct "+str(min_struct) \
                    + " --align_min_bright "+str(min_bright) \
                    + " --noise "+str(noise)+" --drizzle "+str(drizzle)
            #command =  "PlanetarySystemStacker "
            #command += temp_dir.name
            #command += " --out_format fits -b 4 -s "+str(int(stack_percent))
            #command += " -a 52 -w 20 --normalize_bright"

            #   Run Planetary System Stacker command
            # subprocess.run([command], shell=True, text=True)
            subprocess.run([command], shell=True)

            #   Create a link to the stack images in a temporary directory
            os.symlink(
                temp_dir.name+'_pss.tiff',
                #temp_dir.name+'_pss.'+filetype,
                temp_dir_stack.name+'/'+str(j).rjust(digits, '0')+"_stack."+filetype
                )
            j += 1

            #   Create new temporary directory
            temp_dir = tempfile.TemporaryDirectory()

    #   Create new file collection
    if bfits:
        #ifc = ccdp.ImageFileCollection(filenames=path_stacked_imgs)
        ifc = ccdp.ImageFileCollection(temp_dir_stack.name)

        #   Number of files
        nfiles = len(ifc.files)
    else:
        files, nfiles = aux.mkfilelist(
            temp_dir_stack.name,
            formats=formats,
            addpath=True,
            sort=True,
            )

        ifc = imread_collection(files)


    #   Scale trimming indices to account for drizzling
    if drizzle != 'Off':
        drizzle = float(drizzle[:-1])
        ys_cut = int(ys_cut*drizzle)
        ye_cut = int(ye_cut*drizzle)
        xs_cut = int(xs_cut*drizzle)
        xe_cut = int(xe_cut*drizzle)


###
#   Get the best images, using the Planetary System Stacker
#
if best_img and not stack:
    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    #   Load configuration for the Planetary System Stacker
    configuration = Configuration()
    configuration.initialize_configuration(read_from_file=False)

    #   Prepare a list for images in ''step''
    path_list = []

    for i, img_path in enumerate(files):
        #   Add image to the list
        path_list.append(img_path)

        #   If set is reach, find best images
        if ((i != 0 and i%step == 0) or
            (i+1 == len(files) and i%step != 0 and add_last_best)):
            #   Get images as a frames collection
            frames = Frames(configuration, path_list, type='image')

            #   Rank images
            ranked_frames = RankFrames(frames, configuration)
            ranked_frames.frame_score()

            #   Identify best images
            best_indices, _, _ = ranked_frames.find_best_frames(nimgs, window)

            #print ("\nIndices of the best frames "+str(best_indices)+" in step number "
                #+str(int(i/step)))
            for indice in best_indices:
                #   Path to best image
                path_best = path_list[indice]

                #   Get image file name
                filename = str(path_best).split('/')[-1]

                #   Create a link to the best images in the temporary directory
                os.symlink(pwd+'/'+path_best, temp_dir.name+'/'+filename)

                #   Reset image list
                path_list = []


    #   Create new file collection
    if bfits:
        ifc = ccdp.ImageFileCollection(temp_dir.name)

        #   Number of files
        nfiles = len(ifc.files)
    else:
        files, nfiles = aux.mkfilelist(
            temp_dir.name,
            formats=formats,
            addpath=True,
            sort=True,
            )

        ifc = imread_collection(files)


###
#   Cut pictures and postprocess them
#
if bfits:
    #   Loop over and trim all images
    i = 0
    #img_list = []
    for img_ccd, fname in ifc.ccds(ccd_kwargs={'unit': 'adu'}, return_fname=True):
        #   Write status to console
        sys.stdout.write("\rApply cuts to images %i/%i" % (i+1, nfiles))
        sys.stdout.flush()
        i += 1

        #   Trim images
        img_out = ccdp.trim_image(img_ccd[ys_cut:-ye_cut, xs_cut:-xe_cut])

        #   Postprocess image
        img_out.data = postprocess.postprocess(
            img_out.data,
            bit_depth,
            global_histo_equalization,
            local_histo_equalization,
            disk_size,
            adaptive_histo_equalization,
            clip_limit,
            log_cont_adjust,
            log_gain,
            gamma_adjust,
            gamma,
            contrast_stretching,
            upper_percentile,
            lower_percentile,
            postprocess_layers,
            )


        ###
        #   Prepare array with RGB image
        #
        if mkRGB and out_format not in [".FIT",".fit",".FITS",".fits"]:
            #   Get shape of the trimmed image
            out_shape = img_out.data.shape

            #   Prepare array
            rgb_img = np.zeros(
                (out_shape[0], out_shape[1], 4),
                dtype='uint8',
                )

            #   Scale data, convert data to 8bit range and add data to the array
            rgb_img[:,:,0] = img_out.data * r_scale / 2**(bit_depth) * 255
            rgb_img[:,:,1] = img_out.data * g_scale / 2**(bit_depth) * 255
            rgb_img[:,:,2] = img_out.data * b_scale / 2**(bit_depth) * 255
            rgb_img[:,:,3] = 255


        ###
        #   Write postprocess image
        #
        new_name = os.path.basename(fname).split('.')[0]
        if out_format in [".tiff", ".TIFF"]:
            if mkRGB:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    rgb_img,
                    )
            else:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    img_out.data,
                    )
        elif out_format in [".jpg", ".jpeg", ".JPG", ".JPEG"]:
            if mkRGB:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    rgb_img[:,:,0:3],
                    )
            else:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    img_out.data,
                    )
        elif out_format in [".FIT",".fit",".FITS",".fits"]:
            img_out.write(post_path / fname, overwrite=True)
        else:
            print('Error: Output format not known :-(')


        #   Image normalizatiom by means of ImageMagick
        if norm_image_magick:
            command = ('convert -normalize '+os.path.join(path_out, 'postpro', new_name)+out_format+' '+os.path.join(path_out, 'postpro', new_name)+out_format)
            #print(command)
            subprocess.run(
            [command],
            shell=True,
            #text=True,
            #capture_output=True,
            )
else:
    for i, fname in enumerate(files):
        #   Write status to console
        id_c = i+1
        sys.stdout.write("\rApply shift to image %i/%i" % (id_c, nfiles))
        sys.stdout.flush()

        #   Trim images
        if xe_cut != 0 and ye_cut != 0:
            img_i = ifc[i][ys_cut:-ye_cut, xs_cut:-xe_cut]
        elif xe_cut == 0 and ye_cut != 0:
            img_i = ifc[i][ys_cut:-ye_cut, xs_cut:]
        elif xe_cut != 0 and ye_cut == 0:
            img_i = ifc[i][ys_cut:, xs_cut:-xe_cut]
        else:
            img_i = ifc[i][ys_cut:, xs_cut:]

        #   Postprocess image
        img_i = postprocess.postprocess(
            img_i,
            bit_depth,
            global_histo_equalization,
            local_histo_equalization,
            disk_size,
            adaptive_histo_equalization,
            clip_limit,
            log_cont_adjust,
            log_gain,
            gamma_adjust,
            gamma,
            contrast_stretching,
            upper_percentile,
            lower_percentile,
            postprocess_layers,
            )

        ###
        #   Prepare array with RGB image
        #
        if mkRGB and out_format not in [".FIT",".fit",".FITS",".fits"]:
            #   Get shape of the trimmed image
            out_shape = img_i.shape

            #   Prepare array
            rgb_img = np.zeros(
                (out_shape[0], out_shape[1], 4),
                dtype='uint8',
                )

            #   Scale data, convert data to 8bit range and add data to the array
            rgb_img[:,:,0] = img_i * r_scale / 2**(bit_depth) * 255
            rgb_img[:,:,1] = img_i * g_scale / 2**(bit_depth) * 255
            rgb_img[:,:,2] = img_i * b_scale / 2**(bit_depth) * 255
            rgb_img[:,:,3] = 255



        ###
        #   Write postprocess image
        #
        new_name = os.path.basename(fname).split('.')[0]
        default_formats = [
            ".jpg",
            ".jpeg",
            ".JPG",
            ".JPEG",
            ".FIT",
            ".fit",
            ".FITS",
            ".fits",
            ]
        if out_format in [".tiff", ".TIFF"]:
            if mkRGB:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    rgb_img,
                    )
            else:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    img_i,
                    )
        elif out_format in default_formats:
            if mkRGB:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    rgb_img[:,:,0:3],
                    )
            else:
                imsave(
                    os.path.join(path_out, 'postpro', new_name)+out_format,
                    img_i,
                    )
        else:
            print('Error: Output format not known :-(')


        #   Image normalizatiom by means of ImageMagick
        if norm_image_magick:
            command = ('convert -normalize '+os.path.join(path_out, 'postpro', new_name)+out_format+' '+os.path.join(path_out, 'postpro', new_name)+out_format)
            #print(command)
            subprocess.run(
            [command],
            shell=True,
            #text=True,
            #capture_output=True,
            )

sys.stdout.write("\n")

