#!/usr/bin/env python3

import os, pylab, time, sys, matplotlib
from numpy import matrix
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from decimal import Decimal
import argparse

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def poly3(x,m1,m2,m3,b):
    p3 = m1*(x**3) + m2*(x**2) + m3*x + b
    return p3

def gauss(a,x,x0,sigma):
    g = a * np.exp(-(x-x0)**2/(2*sigma**2))
    return g

def gaussian(initial,x):
    A, x1, sig1, cont = initial
    g = A * np.exp(-(x-x1)**2/(2*sig1**2)) + cont
    return g

def gaussian_double(initial,x):
    A,x1,sig1,B,x2,sig2,cont = initial
    g2 = A*np.exp(-(x-x1)**2/2*sig1**2) + B*np.exp(-(x-x2)**2/2*sig2**2) + cont
    return g2

def res_gauss(p,x,y):
    return gaussian(p,x) - y

def res_gauss_double(p,x,y):
    return gaussian_double(p,x) - y

def fit_gaussian_line(wave, flux, line_center, **kwargs):
    from scipy import optimize

    wave = np.array(wave)
    flux = np.array(flux)

    # only consider close pixels for the "continuum"
    hxrange = 50
    mask = np.where((line_center - hxrange < wave) & (wave < line_center + hxrange))
    x = wave [mask]
    y = flux [mask]
    # take median of lowest points (cont. for emission lines)
    y_sorted = np.sort(y)
    cut_thres = 0.5
    cut_index = int(len(y_sorted)*cut_thres)
    y_sorted_cut = y_sorted [:cut_index]
    cont = np.median(y_sorted_cut)

    # first guess
    line_ew = 1.5
    sigma = 100
    p0 = [line_ew, line_center, sigma, cont]

    # fit gaussian
    p1, conv = optimize.leastsq(res_gauss, p0[:], args=(x,y))
    print(p1)

    if "double" in kwargs:
        line_center_1 = line_center - 2
        line_center_2 = line_center + 2
        sigma_1 = sigma
        sigma_2 = sigma
        line_ew_1 = line_ew
        line_ew_2 = line_ew
        p0 = [line_ew_1, line_center_1, sigma_1,
              line_ew_2, line_center_2, sigma_2, cont]
        p1, conv = optimize.leastsq(res_gauss_double,p0[:],args=(x,y))

    if (conv > 3): #or (abs(p1[1] - p0[1]) > 4):
        print("Gauss fit not converged", conv, p1)
    return p1

def click_point(event,linelist_tmp,fig1,lines_found_idx,spec,
                lines_found_int,lines_identified_wave,lines_identified_pixel,sname_plot_select):

    if event.button == 1:
        #print('left click on '+str(event.xdata),str(event.ydata))
        #fig1.suptitle('last click: left')

        if len(linelist_tmp) > 0:
            # compute lclick distance to each identified line
            distance = []
            for i in range(len(lines_found_idx)):
                distance_i = np.sqrt(((event.xdata - lines_found_idx[i]) / \
                    len(spec))**2 + ((event.ydata - lines_found_int[i]) / \
                        ((max(spec) + 0.1*max(spec))*(16/9)))**2)
                distance.append(distance_i)
                #distance.append(np.sqrt((event.xdata)**2+(event.ydata)*ratio)**2)

            # select the closest line
            index_mindistx = distance.index(min(distance))
            lines_identified_wave.append(linelist_tmp[0])
            lines_identified_pixel.append(lines_found_idx[index_mindistx])

            plt.text(lines_found_idx[index_mindistx] + len(spec)/100, lines_found_int[index_mindistx],\
                str(linelist_tmp[0]),va="center",ha="left")
            linelist_tmp.pop(0)

            fig1.suptitle(r'$\lambda$ = ' + str( linelist_tmp[0]) + r' $\AA$',color='r',fontsize=20)
            plt.plot(lines_found_idx[index_mindistx], lines_found_int[index_mindistx], 'bo')

            #plt.plot([lines_found_idx[index_mindistx],event.xdata],[lines_found_int[index_mindistx],event.ydata],'b-')
            event.canvas.draw()

            lines_found_idx.pop(index_mindistx)
            lines_found_int.pop(index_mindistx)
        else:
            fig1.suptitle('All lines done')
            time.sleep(1.0)
            fig1.suptitle("")
            plt.plot(lines_found_idx[index_mindistx],lines_found_int[index_mindistx],'bo')
            plt.savefig(sname_plot_select,bbox_inches="tight")
            plt.close()

    if event.button == 3:
        #print('right click on '+str(event.xdata),str(event.ydata))
        if len(linelist_tmp) > 0:
            linelist_tmp.pop(0)
            #fig1.title('spec_intensity_1d',color='g')
            fig1.suptitle(r'$\lambda$ = '+str( linelist_tmp[0])+ r' $\AA$',color='r',fontsize=20)
            event.canvas.draw()
        else:
            fig1.suptitle('All lines done')
            time.sleep(1.0)
            fig1.suptitle("")
            plt.savefig(sname_plot_select,bbox_inches="tight")
            plt.close()

def press_button(event,fig1,sname_plot_select):
    if event.key == 'q':
        fig1.suptitle("")
        plt.savefig(sname_plot_select,bbox_inches="tight")
        plt.close()

def fit_disperison(lines_identified_pixel,lines_identified_wave,order):
    '''polynomial fit to the wavelength (pixel) relation'''
    if len(lines_identified_pixel) < 4:
        print(bcolors.FAIL+"   [ERROR] Not enough lines selected for a fit. Try again!"+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
        sys.exit()
    print(bcolors.BOLD+"   Fit selected lines ("+str(len(lines_identified_pixel))+") with 3rd order polynomial "+bcolors.ENDC)
    m1,m2,m3,b = pylab.polyfit(lines_identified_pixel,lines_identified_wave,3)
    #print(m1,m2,m3,b)
    print(bcolors.OKGREEN+"   Fit: "+'y = '+str('%.2E' % Decimal(m1))+" x^3$ + "+\
        str('%.2E' % Decimal(m2))+" x^2$ + "+str('%.2E' % Decimal(m3))+"x + " +\
            str('{:7.2f}'.format(b))+bcolors.ENDC)
    return m1,m2,m3,b

def plot_disp(lines_identified_pixel,lines_identified_wave,m1,m2,m3,b,spec_pixel,sname_plot_disp):
    '''plot the dispersion relation'''
    fig2 = plt.figure(num='Wavelength calibration fit')

    plt.ylabel(r'$\lambda\,[\AA]$')
    plt.xlabel('Pixel')

    fig2.suptitle('$\lambda$ = '+str('%.2E' % Decimal(m1))+r" $\times$ pxl$^3$ + "+\
        str('%.2E' % Decimal(m2))+r" $\times$ pxl$^2$ + "+str('%.2E' % Decimal(m3))+\
            r" $\times$ pxl + " + str('{:7.2f}'.format(b)),fontsize=13)
    #fig2.canvas.mpl_connect('key_press_event',press_key)

    plt.plot(spec_pixel,poly3(spec_pixel,m1,m2,m3,b),'r--',linewidth=0.5)
    plt.plot(lines_identified_pixel,lines_identified_wave,'go')

    print(bcolors.BOLD+"   Fit is displayed in external window"+bcolors.ENDC)
    print(bcolors.BOLD+"   This plot will be saved as "+bcolors.OKBLUE+sname_plot_disp+bcolors.ENDC)
    print(bcolors.OKBLUE+"   [ACTION REQUIRED] Close window to continue"+bcolors.ENDC)
    plt.tight_layout()
    plt.savefig(sname_plot_disp)
    plt.show('Wavelength calibration fit')
    plt.clf()
    #plt.close()

def ifilter_median(array,limit):
    from scipy.ndimage import median_filter
    img_filtered = np.array(median_filter(array,size=limit)).astype(np.float32)
    return img_filtered

def find_rspec(iname):
    ''' this only shows a comparison spectrum for now'''
    from matplotlib.colors import LogNorm

    with fits.open(iname) as hdu_list:
        image_data = hdu_list[0].data
        hdu_list.info()

    #image_data = ifilter_median(image_data,100)
    print(image_data.shape)

    #image_data = fits.getdata(iname)
    plt.imshow(image_data, cmap='gray', norm=LogNorm(),
               vmin=5000,vmax=8000)
    #plt.imshow(image_data, cmap='gray')

    plt.colorbar()
    plt.tight_layout()

    #NBINS = 1000
    #histogram = plt.hist(image_data.flatten(), NBINS)

    plt.show()

def est_cont_rough(spec_tmp):
    from scipy.ndimage import gaussian_filter1d

    wsize = 50
    cutmedian = 0.9
    spec_smooth_median = [np.median(spec_tmp[i-wsize:i+wsize]) for i in range(0,len(spec_tmp))]

    spec_smooth_median = np.array(spec_smooth_median)
    mask = np.isnan(spec_smooth_median)
    spec_smooth_median[mask] = np.interp(np.flatnonzero(mask),
                                      np.flatnonzero(~mask),
                                      spec_smooth_median[~mask])
    spec_smooth_median_gauss = gaussian_filter1d(spec_smooth_median, 13)

    return spec_smooth_median_gauss

def find_emission_lines(spec_int,emission_lines_to_find,emissionlinewidth):
    '''try to find emission line maxima'''
    spec_int = list(spec_int)
    spec_x = range(0,len(spec_int))
    lines_found = []
    gauss_inputs = []
    for i in range(0,emission_lines_to_find + 30):
        max_index = spec_int.index(max(spec_int))
        max_value = max(spec_int)

        #kwargs_fitgauss = {"double":"test"}
        #kwargs_fitgauss = {}
        #gauss_input = fit_gaussian_line(spec_x, spec_int, spec_x[max_index],
        #                                **kwargs_fitgauss)
        #gauss_inputs.append(gauss_input)
        gauss_inputs.append(None)

        spec_int = [spec_int[j] - gauss(max_value,j,max_index,emissionlinewidth) for j in range(0,len(spec_int))]
        #spec_int = [spec_int[j] - gaussian(gauss_input,j) for j in range(0,len(spec_int))]

        lines_found.append([max_index,max_value])

    lines_found = np.array(lines_found).T
    #print(lines_found)
    return lines_found, gauss_inputs

def find_emission_lines_v2(spec_tmp,emission_lines_to_find,emissionlinewidth):
    '''try to find emission line maxima'''
    lines_found = []
    winsize = 20
    print("winsize",winsize)
    #lines_strong = [[i,spec_tmp[i]] for i in range(0,len(spec_tmp)) if spec_tmp[i] > np.median(spec_tmp[i-winsize:i+winsize]) ]
    lines_strong = []
    for i in range(0,len(spec_tmp)):
        if spec_tmp[i] > 2*np.median(spec_tmp[i-winsize:i+winsize]):
            lines_strong.append([i, spec_tmp[i]])
    lines_strong = np.array(lines_strong).T

    return lines_strong

def identify_lines(spec_pixel,spec,lines_found,linelist_tmp,sname_plot_select,**kwargs):
    '''plot calibration spectrum and identified lines in pixel coordinates'''

    lines_found_idx = list(lines_found[0])
    lines_found_int = list(lines_found[1])
    fig1 = plt.figure(num='Select line for the given wavelength',
                      figsize=(16, 9))
    axes = plt.gca()

    axes.set_xlabel('Pixel')
    axes.set_ylabel('Counts')
    axes.set_yscale('log')
    axes.plot(spec_pixel,spec,"b-", linewidth=0.3)
    axes.plot(lines_found[0],lines_found[1],"ro",linewidth=5)

    if "gauss" in kwargs:
        for item in kwargs["gauss"]:
            wcen = item[1]
            hxrange = 10
            mask = np.where((wcen - hxrange < spec_pixel) & (spec_pixel < wcen + hxrange))
            spec_pixel_gauss = spec_pixel [mask]
            axes.plot(spec_pixel_gauss, gaussian(item,spec_pixel_gauss), color='green')

    if "dgauss" in kwargs:
        for item in kwargs["dgauss"]:
            wcen = item[1]
            hxrange = 10
            mask = np.where((wcen - hxrange < spec_pixel) & (spec_pixel < wcen + hxrange))
            spec_pixel_gauss = spec_pixel [mask]
            axes.plot(spec_pixel_gauss, gaussian_double(item,spec_pixel_gauss), color='green')

    axes.set_xlim([0,len(spec)])
    #axes.set_ylim([0,max(spec)+0.1*max(spec)])
    print(min(spec))
    spec_range = max(spec) - min(spec)
    axes.set_ylim([0.9*min([abs(i) for i in spec]),max(spec)+0.1*spec_range])

    fig1.suptitle(r'$\lambda$ = '+str( linelist_tmp[0])+ r' $\AA$',color='r',fontsize=20)

    # ratio = (max(lines_found_idx)-min(lines_found_idx)) / (max(lines_found_int)-min(lines_found_int))
    lines_identified_pixel = []
    lines_identified_wave = []

    # manually identify selected emission lines with their wavelength
    fig1.canvas.mpl_connect('button_press_event',
                            lambda event: click_point(event,linelist_tmp,fig1,
                                                      lines_found_idx,spec,
                                                      lines_found_int,lines_identified_wave,lines_identified_pixel,
                                                      sname_plot_select))
    fig1.canvas.mpl_connect('key_press_event',
                            lambda event: press_button(event,fig1,sname_plot_select))

    print(bcolors.BOLD+"   Emissions lines identified, calibration plot opened"+bcolors.ENDC)

    print(bcolors.OKBLUE+"   [ACTION REQUIRED] Select matching emission line with displayed wavelength"+bcolors.ENDC)
    print(bcolors.OKBLUE+"   [ACTION REQUIRED] Left Click: select line"+bcolors.ENDC)
    print(bcolors.OKBLUE+"   [ACTION REQUIRED] Right Click: skip current wavelength"+bcolors.ENDC)
    print(bcolors.OKBLUE+"   [ACTION REQUIRED] Select at least 4 lines"+bcolors.ENDC)
    print(bcolors.OKBLUE+"   [ACTION REQUIRED] Press 'q' or close plot window to finish selection"+bcolors.ENDC)

    plt.show()
    plt.clf()
    plt.close('all')

    return lines_identified_pixel, lines_identified_wave


def main():

    description = '''\
    description:
      Find position of emission lines in calibration spectrum and derive
      dispersion relation by fitting a 3rd order polynomial.
        '''
    epilog = '''\
    examples:
      ./1_findcaliblines.py -arc Hg_1s_II.FIT -rsc 104 114 -rbg 0 20
      ./1_findcaliblines.py -arc Hg_1s_II.FIT -rsc 520 560 -rbg 0 200 -sc star.fit
      ./1_findcaliblines.py -arc Hg_1s_II.FIT -rsc 104 114 -rbg 0 20 -lw 5 -nl 30
        '''

    parser = argparse.ArgumentParser(
          prog='1_findcaliblines.py',
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description=description,
          epilog=epilog)

    parser.add_argument('-arc',"--arclamp",
                        nargs=1, default=['calib_wave.FIT'],
                        help="calibration (arc) image")

    parser.add_argument('-sc',"--science",
                        nargs=1,
                        help="science image to find best extraction region")

    parser.add_argument('-rsc',"--region_science",
                        nargs=2, default=[495,590],type=int,
                        help="region containing the science spectrum")

    parser.add_argument('-rbg',"--region_background",
                        nargs=2, default=[0,200],type=int,
                        help="region containing the background spectrum")

    parser.add_argument('-lw',"--linewidth",
                        nargs=1, default=[1.5],type=float,
                        help="first guess for the emission line width in pixel")

    parser.add_argument('-nl',"--number_lines",
                        nargs=1, default=[35],type=int,
                        help="number of emission lines to search for")

    parser.add_argument('-wr',"--wave_range",
                        nargs=2, default=[3000,10000],type=float,
                        help="maximum expected wavelength range; lines outside are excluded")

    #parser.add_argument('argv', nargs='*')
    args = parser.parse_args()

    # name of the file with the wavelength calibration spectrum
    calibFileName  =   args.arclamp[0]

    # region (rows on the image) containing the calibration spectrum
    specRegionStart = args.region_science[0]
    specRegionEnd   = args.region_science[1]

    if args.science:
        name_science = args.science[0]
        # show science frame so relevant rows can be identified
        find_rspec(name_science)

    # background region (rows on the image), which needs to be outside of the slits
    bgRegionStart   = args.region_background[0]
    bgRegionEnd     = args.region_background[1]

    # width of emission lines in pixel
    emissionlinewidth = args.linewidth[0]

    # names of output files
    sname_spec = "calibration_spectrum.dat"
    sname_plot_disp = "calibration_fit.pdf"
    sname_plot_select = "calibration_selection.pdf"

    # expected lines in calibration (arc) spectrum

    # HgAr (standard fluorescent tube)
    linelist_hgar = [3650.158,4046.56,4077.837,4358.33,5460.75,5769.610,
                5790.670,6965.431,7067.218,7147.042,7272.936,7383.980,
                7503.86,7635.44,7724.20,7948.176,8006.157,8103.693,
                8115.81,8264.522,8408.210,8521.442,8667.944]

    # NeAr calibration lamp (Sheylak Alpy)
    linelist_near = [4072.0043,4103.91,4131.73,4158.59,4164.180,4181.884,
                 4198.317,4200.675,4237.2195,4259.361,4266.286,
                 4272.168,4277.5279,4300.101,4331.25,4333.56,4348.11,
                 4379.74,4426.01,4430.18,4448.8791,4471.586,4474.77,
                 4481.83,4510.733,4545.08,4579.39,4589.93,4609.6,
                 4637.2324,4657.94,4702.316,4726.91,4732.08,4735.93,
                 4764.89,4806.07,4847.9,4879.9,4889.06,4933.2089,
                 4965.12,5017.16,5062.07,5090.55,5141.81,5145.36,
                 5162.284,5187.746,5221.271,5330.7775,5341.0938,
                 5343.2834,5400.5616,5421.352,5439.97,5451.65,
                 5495.872,5506.112,5558.702,5572.465,5606.732,
                 5650.704,5739.520,5748.2985,5764.4188,5820.1558,
                 5834.263,5852.4878,5881.895,5888.584,5902.4623,
                 5912.085,5928.813,5944.834,6043.223,6059.372,
                 6074.3376,6096.163,6114.9232,6143.0627,6163.5937,
                 6266.4952,6304.7893,6334.4276,6382.9914,6402.248,
                 6416.307,6506.5277,6532.8824,6598.9528,6717.043,
                 6752.834]

    # NeXe "party" lamp
    linelist_nexe = [4500.978,4524.681,4582.747,4624.275,4671.226,
                4690.97,4697.021,4734.152,4792.619,4807.019,4829.708,
                4843.293,4916.51,4923.152,5028.279,5400.562,5852.488,
                5881.895,5944.834,5975.534,6029.997,6074.337,6096.163,
                6143.063,6163.59,6182.42,6217.281,6266.495,6304.789,
                6318.062,6334.428,6382.991,6402.248]

    # HeAr lamp used for EFOSC2
    linelist_hear = [3187.743, 3393.73, 3464.14, 3520., 3545.58, 3559.51, 3718.21,
                3729.29, 3737.89, 3780.84, 3819.6074, 3850.57, 3888.646, 3928.62,
                3948.979, 3964.727, 4026.189, 4044.418, 4072.2, 4131.73,
                4158.59, 4200.674, 4259.361, 4277.55, 4300.4, 4333.561,4387.9296,
                4426.01, 4471.477, 4510.733, 4545.08, 4579.39, 4657.94,
                4713.143, 4764.89, 4806.07, 4879.9, 4921.929, 4965.12,
                5015.675, 5047.738, 5187.746, 5221.27, 5495.872, 5572.548, 5606.732,
                5650.703, 5875.618, 6114.923, 6172.278, 6416.3, 6678.149,
                6752.832, 6871.29, 6965.43, 7065.188, 7107.496, 7125.8,
                7147.041, 7206.986, 7272.936, 7281.349, 7311.71, 7353.316,
                7372.118, 7383.98, 7503.867, 7514.651, 7635.105, 7670.04,
                7723.8, 7891.075, 7948.175, 8006.156, 8012., 8014.786,
                8053.307, 8103.692, 8110., 8115.311, 8264.521, 8408.21,
                8424.647, 8521.441, 8605.78, 8620.47, 8667.943, 8761.72,
                9075.42, 9122.966, 9194.68, 9224.498, 9291.58, 9354.218,
                9657.784, 9784.501, 10830.3]

    # missing lines: 3819.6074, 4387.9296, 5047.738
    # standard arc lamp: HgAr fluorescent tube
    linelist = linelist_hgar

    if "NEAR" in calibFileName.upper():
        linelist = linelist_near
    if "NEXE" in calibFileName.upper():
        linelist = linelist_nexe
    if "HEAR" in calibFileName.upper():
        linelist = linelist_hear

    # limit to lines in expected range
    linelist = [i for i in linelist if ((i>min(args.wave_range))&(i<max(args.wave_range))) ]
    #print("len(linelist)",len(linelist))

    # number of emission lines to search for
    if args.number_lines[0] == 'nlines':
        emission_lines_to_find = len(linelist)
    else:
        emission_lines_to_find = args.number_lines[0]

    matplotlib.rcParams['pdf.fonttype'] = 42

    print(bcolors.BOLD+"   Check input data"+bcolors.ENDC)

    if os.path.isfile(calibFileName)==False:
        print(bcolors.FAIL+"   [ERROR] Calibration files doesn't exist."+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
        sys.exit()

    print(bcolors.BOLD+"   Read fits file"+bcolors.ENDC)

    calibFile     = fits.open(calibFileName)
    calibFileData = calibFile[0].data

    # get dimensions of the data
    calibFileDataColumns = len(calibFileData[0]) # number of columns in the file
    calibFileDataLines   = len(calibFileData)    # number of lines in the file
    N = calibFileDataColumns

    if bgRegionStart > calibFileDataLines or bgRegionEnd > calibFileDataLines:
        print(bcolors.FAIL+"   [ERROR] Background region outside of image."+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
        sys.exit()

    if specRegionStart > calibFileDataLines or specRegionEnd > calibFileDataLines:
        print(bcolors.FAIL+"   [ERROR] Background region outside of image."+bcolors.ENDC)
        print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
        sys.exit()

    print(bcolors.BOLD+"   Extract spectral data"+bcolors.ENDC)

    # take average over extraction region set by -rsc
    spec = []
    for i in range(0,N):
        spec_tmp, bg_tmp = 0,0
        for j in range(specRegionStart,specRegionEnd):
            spec_tmp += calibFileData[j][i]
        for k in range(bgRegionStart,bgRegionEnd):
            bg_tmp += calibFileData[k][i]
        spec_tmp /= (specRegionEnd-specRegionStart)
        bg_tmp /= (bgRegionEnd-bgRegionStart)
        #spec.append(spec_tmp-bg_tmp)
        spec.append(spec_tmp)

    spec_pixel  = np.arange(len(spec))
    spec_intensity_1d = np.asarray(spec)

    # try to identify emission line maxima in pixel coordinates
    print(bcolors.BOLD+"   Search for "+str(emission_lines_to_find)+\
          " emissions lines"+bcolors.ENDC)
    spec_tmp = spec.copy()
    spec_tmp_smooth = est_cont_rough(spec_tmp)
    spec_tmp_s = spec_tmp / spec_tmp_smooth

    lines_found, gauss_inputs = find_emission_lines(spec_tmp_s,
                                                    emission_lines_to_find,
                                                    emissionlinewidth)
    #lines_found = [lines_found[i] for i in range(0,len(lines_found)) if gauss_inputs[i][1] != None]
    #gauss_inputs = [gauss_inputs[i] for i in range(0,len(lines_found)) if gauss_inputs[i][1] != None]

    #lines_found = find_emission_lines_v2(spec_tmp,emission_lines_to_find,emissionlinewidth)

    # plot calibration spectrum and identified lines in pixel coordinates
    linelist_tmp = linelist

    #kwargs_id = {"gauss":gauss_inputs}
    #kwargs_id = {"dgauss":gauss_inputs}
    kwargs_id = {}

    lines_identified_pixel, lines_identified_wave = \
        identify_lines(spec_pixel,spec_tmp_s,
                      lines_found,linelist_tmp,
                      sname_plot_select,**kwargs_id)

    m1, m2, m3, b = fit_disperison(lines_identified_pixel,lines_identified_wave,3)
    plot_disp(lines_identified_pixel,lines_identified_wave,
              m1,m2,m3,b,spec_pixel,sname_plot_disp)

    print(bcolors.BOLD+"   Write calibrated calibration spectrum to "+\
        bcolors.OKBLUE+sname_spec+bcolors.ENDC)

    os.system('touch '+sname_spec)
    with open(sname_spec,'w') as specout:
        for i in range(len(spec_pixel)):
            specout.write(str(poly3(spec_pixel[i],m1,m2,m3,b))+'  '+str(spec_intensity_1d[i])+"\n")
        print(bcolors.OKGREEN+"   DONE"+bcolors.ENDC)

    os.system('touch '+sname_spec+'.norm')
    with open(sname_spec+'.norm','w') as specout:
        for i in range(len(spec_pixel)):
            specout.write(str(poly3(spec_pixel[i],m1,m2,m3,b))+'  '+str(spec_tmp_s[i])+"\n")
        print(bcolors.OKGREEN+"   DONE"+bcolors.ENDC)

if __name__ == "__main__":
    main()
