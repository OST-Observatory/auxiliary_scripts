#!/usr/bin/env python3

import os, time, sys, matplotlib
from astropy.io import fits
import numpy as np
from scipy.interpolate import interp1d
import argparse
import matplotlib.pyplot as plt

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def plot_spectra_reduced(spec_science,spec_flat,
                         wave_grid,plotident,savename_sci_pdf):

    fig1 = plt.figure(figsize=(20,10))
    axes = plt.gca()

    font = {'family': 'serif',
            'color':  'red',
            'weight': 'normal',
            'size': 15}

    # setting plot labels 
    axes.set_xlabel(r'$\lambda\,[\AA]$')
    axes.set_ylabel('Relative flux')

    # setting plot ranges

    ymin_plt  =  min(spec_science[100:-100])
    ymax_plt  = max(spec_science[100:-100])
    yrange_pad = (ymax_plt - ymin_plt)*0.02
    axes.set_ylim(ymin_plt-yrange_pad,
             ymax_plt+yrange_pad)

    xrange_pad = (max(wave_grid) - min(wave_grid)) * 0.01
    axes.set_xlim(min(wave_grid)-xrange_pad,
             max(wave_grid)+xrange_pad)

    # plot the actual data
    axes.plot(wave_grid, spec_science, 'b-')

    # setting plotpositions for ident lines
    if plotident == 'yes':
        plotminimum    = ymin_plt - yrange_pad
        plotmaximum    = ymax_plt + yrange_pad
        plotheight     = plotmaximum - plotminimum
        plotmiddleplot = (plotminimum + plotmaximum)/2.0
        plotupperplot  = plotminimum + 0.80*plotheight
        plotlowerplot  = plotminimum + 0.20*plotheight
        plotuppercut1  = plotminimum + 0.70*plotheight
        plotlowercut1  = plotminimum + 0.30*plotheight
        plotuppercut2  = plotminimum + 0.68*plotheight
        plotlowercut2  = plotminimum + 0.32*plotheight

        # interpolate on data to find point for ident
        f2 = interp1d(wave_grid, spec_science)

        lineFile = 'lines.txt'
        lines = open(lineFile,"r")

        # plot idents to figure
        for line in lines:
            frow = line.split()
            if len(frow) == 1:
                print(bcolors.WARNING+"     [WARNING] Broken identification found as '"+line+"', must consist of [wavelenght(s) + name]. Skipping."+bcolors.ENDC)
                continue
            try:
                float(frow[0])
            except ValueError:
                print(bcolors.WARNING+"     [WARNING] Broken identification found as '"+line+"', first entry not a number. Skipping."+bcolors.ENDC)
                continue

            if len(frow) == 2:
                ident_line = frow
                # single ident plot
                if float(ident_line[0]) >= min(wave_grid) and float(ident_line[0]) <= max(wave_grid):
                    if f2(ident_line[0]) <= plotmiddleplot:
                        axes.plot([ident_line[0], ident_line[0]],
                                  [f2(ident_line[0]), plotupperplot],
                                  color='r', linestyle='-', linewidth=1.5)
                        axes.text(ident_line[0], plotupperplot, ident_line[1],
                                  rotation=90, ha='center', va='bottom', fontdict=font)
                    else:
                        axes.plot([ident_line[0], ident_line[0]], [f2(ident_line[0]), plotlowerplot],
                                  color='r', linestyle='-', linewidth=1.5)
                        axes.text(ident_line[0], plotlowerplot, ident_line[1],
                                  rotation=90, ha='center', va='top', fontdict=font)
            # multi ident plot
            if len(frow) > 2:
                points = []
                for i in frow[:-1]:
                    points.append(float(i) )
                pointcenter = sum(points) / float(len(points))
                ident_name = str(frow[-1:])
                # print(points," give ",pointcenter," bei ",frow[-1:])
                if max(points) <= max(wave_grid) and min(points) >= min(wave_grid):
                    if f2(pointcenter) <= plotmiddleplot:
                        axes.plot([pointcenter,pointcenter],[plotupperplot,plotuppercut1],
                                  color='r', linestyle='-', linewidth=1.5)
                        axes.text(pointcenter, plotupperplot, ident_name[2:-2],
                                  rotation=90, ha='center', va='bottom', fontdict=font)
                        for element in points:
                            axes.plot([element,element],[f2(element), plotuppercut2],
                                      color='r', linestyle='-', linewidth=1.5)
                            axes.plot([pointcenter,element],[plotuppercut1, plotuppercut2],
                                       color='r', linestyle='-', linewidth=1.5)
                    if f2(pointcenter) > plotmiddleplot:
                        axes.plot([pointcenter,pointcenter],[plotlowerplot,plotlowercut1],
                                  color='r', linestyle='-', linewidth=1.5)
                        axes.text(pointcenter, plotlowerplot, ident_name[2:-2],
                                  rotation=90, ha='center', va='top', fontdict=font)
                        for element in points:
                            axes.plot([element,element],[f2(element), plotlowercut2],
                                      color='r', linestyle='-', linewidth=1.5)
                            axes.plot([pointcenter,element],[plotlowercut1, plotlowercut2],
                                       color='r', linestyle='-', linewidth=1.5)
        lines.close()

    plt.savefig(savename_sci_pdf,bbox_inches='tight')
    plt.clf()

    # plot the flatfield
    fig2 = plt.figure(figsize=(20,10))

    plt.xlabel(r'$\lambda\,[\AA]$')
    plt.ylabel('Relative flux')

    yoffset = (max(spec_flat)-min(spec_flat)) * 0.05
    plt.ylim([min(spec_flat)-yoffset,max(spec_flat)+yoffset])

    xrange_pad = (float(max(wave_grid)) - float(min(wave_grid))) * 0.01
    plt.xlim(min(wave_grid)-xrange_pad,max(wave_grid)+xrange_pad)

    plt.plot(wave_grid,spec_flat,'r-')

    print(bcolors.BOLD+"   Create flatfield plot "+\
        bcolors.OKBLUE+'spectrum_flat.pdf'+bcolors.ENDC)

    plt.savefig('spectrum_flat.pdf',bbox_inches='tight')

def save_arrays_as_ascii(savename,array_list):
    array_list = [np.array(i).astype(float) for i in array_list]
    spec_save = np.transpose([array_list[0],array_list[1]])
    with open(savename,"w") as savefile:
        np.savetxt(savename,spec_save)
        print(bcolors.BOLD+"   Saved to "+\
            bcolors.OKBLUE+savename+bcolors.ENDC)

def errormessage(text):
    print(bcolors.FAIL+"   [ERROR] "+text+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit()

def average_rows(image,row_start,row_end):
    '''sum rows in science extraction region to produce 1d-spectrum'''
    spec_1d = []
    # i = column number (x-axis)
    for i in range(0,len(image[0])):
        # j = row number (y-axis)
        intensities_col = [image[row][i] for row in range(row_start,row_end)]
        intensity_average = np.mean(intensities_col)
        spec_1d.append(intensity_average)
    return spec_1d

def median_rows(image,row_start,row_end):
    '''sum rows in science extraction region to produce 1d-spectrum'''
    spec_1d = []
    # i = column number (x-axis)
    for i in range(0,len(image[0])):
        # j = row number (y-axis)
        intensities_col = [image[row][i] for row in range(row_start,row_end)]
        intensity_median = np.median(intensities_col)
        spec_1d.append(intensity_median)
    return spec_1d

def main():

    ####################################################################
    ############################ Script Parameters #####################
    ####################################################################


    matplotlib.rcParams['pdf.fonttype'] = 42

    description = '''\
    description:
      Reduction of the stellar spectrum. First, the dark frame is 
      subtracted from the spectrum, then the spectrum is divided by 
      the flatfield, and, lastly, the wavelength calibration is 
      performed. There is also the possibility to mark spectral lines 
      in the spectrum.
        '''
    epilog = '''\
    examples:
      ./2_extractspectrum.py -sc star.FIT -df ../darkframes/300s/ -ff ../flats/ -fd ../darkframes/1s/
      ./2_extractspectrum.py -sc star.FIT -df ../darkframes/300s/ -ff ../flats/ -fd ../darkframes/1s/ -rsc 495 590 -rsk 480 490
        '''

    parser = argparse.ArgumentParser(
          prog='2_extractspectrum.py',
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description=description,
          epilog=epilog)

    parser.add_argument('-sc',"--science",
                        nargs=1, default=['star.FIT'],
                        help="file with science frame")
    parser.add_argument('-df',"--dark_dir",
                        nargs=1, default=['../dark/300s/'],
                        help="darkframe directory (for the science frame)")
    parser.add_argument('-ff',"--flatfield_dir",
                        nargs=1, default=['../dark/300s/'],
                        help="flat frame directory")
    parser.add_argument('-fd',"--flatdark_dir",
                        nargs=1, default=['flats/'],
                        help="darkframe directory (for the flat frame)")
    parser.add_argument('-rsc',"--region_science",
                        nargs=2, default=[495,590],type=int,
                        help="region containing the science spectrum")
    parser.add_argument('-rsk',"--region_sky",
                        nargs=4, default=[0,230,800,1000],type=int,
                        help="two sky background regions (inside the slit)\
                        otherwise: background removal")
    parser.add_argument('-wr',"--wave_range",
                        nargs=2, default=['?','?'],
                        help="plot range, set the variables to '?' for automatic resizing")
    parser.add_argument('-pi',"--plot_ident",
                        nargs=1, default=['no'],
                        help="plot idents: ['yes','no']")
    parser.add_argument('-lf',"--line_file",
                        nargs=1, default=['absorption_lines'],
                        help="file containing line identifications")
    parser.add_argument('-m',"--mode",
                        nargs=1, default=['median'],
                        help="image reduction mode: ['mean','median']")

    args = parser.parse_args()

    # file with science spectrum
    science  =   args.science[0]

    # directory of the darkframe for the science spectrum
    dark_dir  =   args.dark_dir[0]

    # flatfield directory
    flatfield_dir  =   args.flatfield_dir[0]

    # directory of the darkframe for the flats
    flatdark_dir  =   args.flatdark_dir[0]

    # region containing the science spectrum
    rows_sci_start = args.region_science[0]
    rows_sci_end   = args.region_science[1]

    # sky background region (inside the slit)
    rows_sky_hi_start   = args.region_sky[0]
    rows_sky_hi_end     = args.region_sky[1]
    rows_sky_lo_start   = args.region_sky[2]
    rows_sky_lo_end     = args.region_sky[3]

    # plot range, set the variables to '?' for automatic resizing
    lambdamin, lambdamax = args.wave_range[0], args.wave_range[1]

    # line identifications ('yes' or 'no')
    plotident = args.plot_ident[0]

    # file containing line identifications
    lineFile   = args.line_file[0]

    ####################################################################
    #### The following parameters usually do not need to be adjusted ###
    ####################################################################

    mode = args.mode[0]

    # file that contains calibrated calibration spectrum
    calib = 'calibration_spectrum.dat'

    # data output names

    savename_sci_pdf = 'star_spectrum.pdf'
    savename_sci_ascii = 'star_spectrum.dat'

    masterdark_output_name = 'master_dark.fit'
    masterflat_output_name = 'master_flat.fit'

    ####################################################################
    ############################     checks      #######################
    ####################################################################

    print(bcolors.BOLD+"   Check input data"+bcolors.ENDC)

    if os.path.isfile(science) == False:
        errormessage("File containing spectrum doesn't exist.")
    if os.path.exists(dark_dir) == False:
        errormessage("Darkframe directory doesn't exist.")
    if os.path.isfile(dark_dir) == True:
        errormessage("Darkframes are linked to file, give the directory instead!")
    if os.path.exists(flatfield_dir) == False:
        errormessage("Flatfield directory doesn't exist.")
    if os.path.isfile(flatfield_dir) == True:
        errormessage("Flatfields are linked to file, give the directory instead!")
    if os.path.exists(flatdark_dir) == False:
        errormessage("Flatdark directory doesn't exist.")
    if os.path.isfile(flatdark_dir) == True:
        errormessage("Flatdarks are linked to file, give the directory instead!")
    if os.path.isfile(calib) == False:
        errormessage("Calibration file doesn't exist.")

    flatdark_list  = os.listdir(flatdark_dir)
    flatfield_list = os.listdir(flatfield_dir)
    darkframe_list = os.listdir(dark_dir)

    if len(flatdark_list) == 0:
        errormessage("Flatdark directory is empty.")
    if len(flatfield_list) == 0:
        errormessage("Flatfield directory is empty.")
    if len(darkframe_list) == 0:
        errormessage("Darkframe directory is empty.")

    with fits.open(science) as science_image:
        science_image_data = science_image[0].data
        size_science_x = len(science_image_data)
        size_science_y = len(science_image_data[1])

    with fits.open(dark_dir+"/"+darkframe_list[0]) as setup_dark_image:
        setup_dark_image_data = setup_dark_image[0].data
        size_dark_x = len(setup_dark_image_data)
        size_dark_y = len(setup_dark_image_data[1])

    with fits.open(flatfield_dir+"/"+flatfield_list[0]) as setup_flat_image:
        setup_flat_image_data = setup_flat_image[0].data
        size_flat_x = len(setup_flat_image_data)
        size_flat_y = len(setup_flat_image_data[1])

    with fits.open(flatdark_dir+"/"+flatdark_list[0]) as setup_flatdark_image:
        setup_flatdark_image_data = setup_flatdark_image[0].data
        size_flatdark_x = len(setup_flatdark_image_data)
        size_flatdark_y = len(setup_flatdark_image_data[1])

    if size_flatdark_x != size_flat_x or size_flatdark_y != size_flat_y:
        print(bcolors.FAIL+"   [ERROR] Flatdarks and Flatfields don't have the same size, check!"+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ("+size_flatdark_x+"x"+size_flatdark_y+") and ("+size_flat_x+"x"+size_flat_y+")"+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
        sys.exit()

    if size_dark_x != size_flat_x or size_dark_y != size_flat_y:
        print(bcolors.FAIL+"   [ERROR] Flatdarks and Flatfields don't have the same size, check!"+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ("+size_flatdark_x+"x"+size_flatdark_y+") and ("+size_flat_x+"x"+size_flat_y+")"+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
        sys.exit()

    if size_dark_x != size_science_x or size_dark_y != size_science_y:
        print(bcolors.FAIL+"   [ERROR] Flatdarks and Flatfields don't have the same size, check!"+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ("+size_flatdark_x+"x"+size_flatdark_y+") and ("+size_flat_x+"x"+size_flat_y+")"+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
        sys.exit()

    ####################################################################
    ####################     data reduction    #########################
    ####################################################################

    print(bcolors.BOLD+"   Average dark (science) from "+str(len(darkframe_list))+" files"+bcolors.ENDC)
    dark_iarray = np.zeros((size_dark_x,size_dark_y,len(darkframe_list)))
    for element in darkframe_list:
        with fits.open(dark_dir+"/"+element) as dark_image:
            dark_image_data = dark_image[0].data
            darkHeader = dark_image[0].header
            dark_iarray[...,darkframe_list.index(element)] = dark_image_data
    if mode == "median":
        dark_iarray_av = np.median(dark_iarray.astype(int), axis=2)
    if mode == "mean":
        dark_iarray_av = np.mean(dark_iarray.astype(int), axis=2)

    hdu = fits.PrimaryHDU(dark_iarray_av,darkHeader)
    hdu.scale('int16', '', bzero=0)
    hduList = fits.HDUList([hdu])
    hduList.writeto(masterdark_output_name, output_verify='exception', overwrite=True)

    print(bcolors.BOLD+"   Average dark (flat) from "+str(len(flatdark_list))+" files"+bcolors.ENDC)
    flatdark_iarray = np.zeros((size_flat_x,size_flat_y,len(flatdark_list)))
    for element in flatdark_list:
        with fits.open(flatdark_dir+"/"+element) as flatdark_image:
            flatdark_image_data = flatdark_image[0].data
            flatdark_iarray[...,flatdark_list.index(element)] = flatdark_image_data
    if mode == "mean":
        flatdark_iarray_av = np.mean(flatdark_iarray.astype(int), axis=2)
    if mode == "median":
        flatdark_iarray_av = np.median(flatdark_iarray.astype(int), axis=2)

    # darkframe correction for flatfields
    print(bcolors.BOLD+"   Average flat from "+str(len(flatfield_list))+" files"+bcolors.ENDC)
    flatfield_iarray = np.zeros((size_flat_x,size_flat_y,len(flatfield_list)))
    for element in flatfield_list:
        with fits.open(flatfield_dir+"/"+element) as flatfield_image:
            flatfield_image_data = flatfield_image[0].data
            flatHeader = flatfield_image[0].header
            flatfield_iarray[...,flatfield_list.index(element)] = \
                flatfield_image_data - flatdark_iarray_av
    if mode == "mean":
        flatfield_iarray_d_av = np.mean(flatfield_iarray.astype(int), axis=2)
    if mode == "median":
        flatfield_iarray_d_av = np.median(flatfield_iarray.astype(int), axis=2)

    hdu = fits.PrimaryHDU(flatfield_iarray_d_av,flatHeader)
    hdu.scale('int16', '', bzero=0)
    hduList = fits.HDUList([hdu])
    hduList.writeto(masterflat_output_name, overwrite=True)

    print(bcolors.BOLD+"   Apply Darkframe correction to science spectrum"+bcolors.ENDC)
    science_iarray = np.asarray(science_image_data)
    science_iarray_d = np.zeros((size_flat_x,size_flat_y,len(flatdark_list)))
    science_iarray_d = science_iarray - dark_iarray_av

    print(bcolors.BOLD+"   Extract spectrum from given row range"+bcolors.ENDC)

    # sum rows in science extraction region to produce 1d-spectrum

    # star spectrum (intensities)
    spec_science = average_rows(science_iarray_d,rows_sci_start,rows_sci_end)

    # flatfield spectrum (intensities)
    spec_flat = average_rows(flatfield_iarray_d_av,rows_sci_start,rows_sci_end)

    # sky background spectrum (intensities)
    spec_sky_hi = median_rows(science_iarray_d,
                               rows_sky_hi_start,rows_sky_hi_end)
    spec_sky_lo = median_rows(science_iarray_d,
                               rows_sky_lo_start,rows_sky_lo_end)

    spec_sky = np.average([spec_sky_hi,spec_sky_lo],
                          axis=0,
                          weights=[abs(rows_sky_hi_start-rows_sky_hi_end),
                                   abs(rows_sky_lo_start-rows_sky_lo_end)])

    print(bcolors.BOLD+"   Apply Flatfield correction to science spectrum"+bcolors.ENDC)

    # normalize the flatfield
    spec_flat = np.asarray(spec_flat) / max(spec_flat)
    # apply flatfield to science spectrum and sky spectrum
    spec_science_f = spec_science / spec_flat
    spec_sky_f = spec_sky / spec_flat
    spec_science_sky = spec_science - np.array(spec_sky)
    spec_science_sky_f = spec_science_f - spec_sky_f

    print(bcolors.BOLD+"   Apply wavelength calibration"+bcolors.ENDC)

    with open(calib,'r') as calibfile:
        wave_grid = []
        for line in calibfile:
            frow = line.split()
            if len(frow) == 0:
                continue
            if lambdamin != '?' and float(frow[0]) < lambdamin:
                spec_science_f = spec_science_f[1:]
                spec_flat = spec_flat[1:]
                continue
            if lambdamax != '?' and float(frow[0]) > lambdamax:
                spec_science_f = spec_science_f[:-1]
                spec_flat = spec_flat[:-1]
                continue
            wave_grid.append(float(frow[0]))

    print(bcolors.BOLD+"   Create spectral plot "+bcolors.OKBLUE+savename_sci_pdf+bcolors.ENDC)
    plot_spectra_reduced(spec_science_sky_f,spec_flat,
                         wave_grid,plotident,savename_sci_pdf)

    save_arrays_as_ascii("spectrum_science_sky_f.dat",[wave_grid,spec_science_sky_f])
    save_arrays_as_ascii("spectrum_science_sky.dat",[wave_grid,spec_science_sky])
    save_arrays_as_ascii("spectrum_science_f.dat",[wave_grid,spec_science_f])
    save_arrays_as_ascii("spectrum_science.dat",[wave_grid,spec_science])
    save_arrays_as_ascii("spectrum_flat.dat",[wave_grid,spec_flat])
    save_arrays_as_ascii("spectrum_sky_f.dat",[wave_grid,spec_sky_f])
    save_arrays_as_ascii("spectrum_sky.dat",[wave_grid,spec_sky])

    print(bcolors.OKGREEN+"   Done"+bcolors.ENDC)

if __name__ == "__main__":
    main()
