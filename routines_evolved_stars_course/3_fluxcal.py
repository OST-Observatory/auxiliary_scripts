#!/usr/bin/env python3

# Credit: http://python4esac.github.io/plotting/specnorm.html
# modified by M.D.

import argparse
import textwrap
import os
import pylab as plt
import numpy as np
from PyAstronomy import pyasl
from scipy.interpolate import splrep,splev
from scipy import interpolate

def shift_to_obs(waves,vrad):
    c = 299792.458 # km/s
    fac_toobs = np.sqrt((1+(vrad/c))/(1-(vrad/c)))
    shifted = np.multiply(waves,fac_toobs)
    return shifted

def onclick(event,wave,flux,**args):
    # when none of the toolbar buttons is activated and the user clicks in the
    # plot somewhere, compute the median value of the spectrum in a 10angstrom
    # window around the x-coordinate of the clicked point. The y coordinate
    # of the clicked point is not important. Make sure the continuum points
    # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
    # if flag 'fix' is set just take the x,y-coordintes of the clicked point
    args = argparse.Namespace(**args)
    toolbar = plt.get_current_fig_manager().toolbar
    if args.window:
        window = float(args.window)
    else:
        window = 1
    if event.button==1 and toolbar.mode=='':
        if args.fix:
            y_event = event.ydata
        else:
            window = ((event.xdata-window/2)<=wave) & (wave<=(event.xdata+window/2))
            y_event = np.median(flux[window])
        plt.plot(event.xdata,y_event,'rs',ms=5,picker=2,label='cont_pnt')
    plt.draw()

def onpick(event,**kwargs):
    # when the user clicks right on a continuum point, remove it
    if event.mouseevent.button==3:
        if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
            event.artist.remove()
            #clicked_points.remove([event.xdata,event.ydata])

def ontype(event,wave,flux,ymin,ymax,yran,filename,args,kwargs={}):
    # when the user hits enter:
    # 1. Cycle through the artists in the current axes. If it is a continuum
    #    point, remember its coordinates. If it is the fitted continuum from the
    #    previous step, remove it
    # 2. sort the continuum-point-array according to the x-values
    # 3. fit a spline and evaluate it in the wavelength points
    # 4. plot the continuum

    args = argparse.Namespace(**args)
    if 'comp' in kwargs:
        wave_comp,flux_comp = kwargs['comp']
    
    global clicked_points
    global cont_pnt_coord

    if event.key=='enter':
        # fit polynomial to clicked points and plot it
        cont_pnt_coord = []
        # store clicked points for easy modification after redo
        clicked_points = []
        
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                cont_pnt_coord.append(artist.get_data())
            elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                artist.remove()
        cont_pnt_coord = np.array(cont_pnt_coord)[...,0]

        sort_array = np.argsort(cont_pnt_coord[:,0])
        x,y = cont_pnt_coord[sort_array].T

        degree = int(args.degree)
        spline = splrep(x,y,k=degree)
        continuum = splev(wave,spline)
    
        plt.plot(wave,continuum,'r-',lw=2,label='continuum')

        clicked_points = list(zip(x,y))

    if event.key=='p':
        # save current clicked points as "clicked_points.dat"
        cont_pnt_coord = []
        clicked_points = []

        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                cont_pnt_coord.append(artist.get_data())
            elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                artist.remove()
        cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
        sort_array = np.argsort(cont_pnt_coord[:,0])
        np.savetxt("clicked_points.dat",cont_pnt_coord[sort_array])
        print("> saved clicked points to 'clicked_points.dat'")

    if event.key=='i':
        # read clicked points from "clicked_points.dat"
        x, y = np.loadtxt("clicked_points.dat",unpack=True,dtype=float)
        cont_pnt_coord = np.array([x,y]).T
        clicked_points = list(zip(x,y))
        print("> read clicked points from 'clicked_points.dat'")
        plt.cla()
        plt.ylim(ymin-0.01*yran,ymax+0.01*yran)
        plt.plot(wave,flux,'k-',
                 lw=0.5,marker='.',ms=0.2)
        for item in clicked_points:
            plt.plot(item[0],item[1],'rs',ms=5,picker=2,label='cont_pnt')

    # when the user hits 'n' and a spline-continuum is fitted,
    # normalise the spectrum
    elif event.key=='n':
        print("> normalising ...")
        global continuum_save
        global continuum_comp_save
        global y_instr
        continuum = None
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                continuum = artist.get_data()[1]
                break
        continuum_save = continuum
        if continuum is not None:
            if args.comparison:
                if args.fluxcal:
                    clicked_wave = np.array(clicked_points).T[0]
                    clicked_flux = np.array(clicked_points).T[1]
                    idx_closest = [min(range(len(wave_comp)), key=lambda i: abs(wave_comp[i]-itm)) for itm in clicked_wave]
                    flux_comp_closet = [flux_comp[i] for i in idx_closest]
                    degree = int(args.degree)
                    print("> fitting spline to comparison ...")
                    spline_comp = splrep(clicked_wave,flux_comp_closet,k=degree)
                    continuum_comp = np.array(splev(wave,spline_comp))
                    continuum_comp_save = continuum_comp
                    #flux_comp_interpol  =  np.interp(wave, wave_comp, flux_comp)
                    print("> dividing continuum by continuum_comp ...")
                    y_instr = continuum/continuum_comp
                    # clip negative values
                    #print("> finding median of instr. profile (to remove negative and small values) ...")
                    #y_instr_median = np.median(y_instr)
                    #y_instr[y_instr<y_instr_median/5e5] = y_instr_median

            plt.cla()
            if not args.fluxcal:
                nymin = min(flux/continuum)
                nymax = max(flux/continuum)
                nyran = nymax - nymin
                plt.ylim(nymin-0.01*nyran,nymax+0.01*nyran)
                plt.plot(wave,flux/continuum,'k-',
                         lw=0.5,marker='.',ms=0.2,label='normalised')
                if args.comparison:
                    #comparison(args.comparison)
                    plt.plot(wave_comp,flux_comp,'k-',
                             lw=0.5,marker='.',ms=0.2,alpha=0.5,
                             label='comparison',c='grey')
                    xmin = min(wave)
                    xmax = max(wave)
                    xlen = xmax-xmin
                    plt.xlim(xmin-0.01*xlen,xmax+0.01*xlen)
            else:
                nymin = min(y_instr)
                nymax = max(y_instr)
                nyran = nymax - nymin
                plt.ylim(nymin-0.01*nyran,nymax+0.01*nyran)
                plt.plot(wave,y_instr,'k-',
                         lw=0.5,marker='.',ms=0.2,label='normalised')

    # when the user hits 'r': clear the axes and plot the original spectrum
    elif event.key=='r':
        plt.cla()
        plt.ylim(ymin-0.01*yran,ymax+0.01*yran)
        plt.plot(wave,flux,'k-',
                 lw=0.5,marker='.',ms=0.2)
        for item in clicked_points:
            plt.plot(item[0],item[1],'rs',ms=5,picker=2,label='cont_pnt')

    # when the user hits 'w': if the normalised spectrum exists, write it to file.
    elif event.key=='w':
        continuum = None
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                continuum = artist.get_data()[1]
                break
        if continuum is not None:
            plt.cla()
            plt.plot(wave,flux/continuum,'k-',
                    lw=0.5,marker='.',ms=0.2,label='normalised')
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                data = np.array(artist.get_data())
                if args.error:
                    errorset = True
                else:
                    errorset = False
                if (errorset == True) & (continuum_save is not None):
                    error_normalized = error/continuum_save
                    data = np.append(data,[error_normalized], axis=0)
                if args.fluxcal:
                    savename = os.path.splitext(filename)[0]+'.instr'
                    savename_fluxcal = os.path.splitext(filename)[0]+'.fluxcal'
                    savedata_fluxcal = np.array([wave,np.array(flux/y_instr)])
                    np.savetxt(savename_fluxcal,savedata_fluxcal.T)
                    print('> saved to',savename_fluxcal)
                    savename_continuum = os.path.splitext(filename)[0]+'.continuum'
                    savedata_continuum = np.array([wave,np.array(continuum_save)])
                    np.savetxt(savename_continuum,savedata_continuum.T)
                    print('> saved to',savename_continuum)
                    savename_continuum_comp = os.path.splitext(filename)[0]+'.continuum_comp'
                    savedata_continuum_comp = np.array([wave,continuum_comp_save])
                    np.savetxt(savename_continuum_comp,savedata_continuum_comp.T)
                    print('> saved to',savename_continuum_comp)
                else:
                    savename = os.path.splitext(filename)[0]+'.nspec'
                np.savetxt(savename,data.T)
                print('> saved to',savename)
                break
    plt.draw()

def convolve(wavec, fluxc, res):
    ''' Apply instr. convolution with constant FWHM \
        FWHM = ((wave[-1]+wave[0])/2) / res'''
    # interpolate to equidistant steps using spline
    if len(wavec)>1:
        stepc   =   abs(wavec[1]-wavec[0])
        if stepc == 0.0:
            stepc =   abs(wavec[2]-wavec[1])
        tck     =   interpolate.splrep(wavec, fluxc, s=0)
        xnew    =   np.arange(wavec[0], wavec[-1], stepc)
        ynew    =   interpolate.splev(xnew, tck, der=0)
        # apply instr. Gauss
        maxsig = 20
        r2, fwhm = pyasl.instrBroadGaussFast(xnew, ynew, res,
            edgeHandling="firstlast", fullout=True, maxsig=maxsig)
        form5   = lambda x: "%.5f" % x
        print("Gaussian kernel FWHM:", form5(fwhm), 'A')
        return xnew, r2
    else:
        print("Warning: Empty range during convolution")
        return wavec, fluxc

def read_spec(filename,**kwargs):
    if 'flat' in kwargs:
        col_flat = int(kwargs['flat'])
        wave,flux=np.loadtxt(filename,skiprows=1,usecols=(0,col_flat),unpack=True,dtype=float)
    else:
        col_wave_obs,col_flux_obs = 0,1
        if 'readcols' in kwargs:
            col_wave_obs, col_flux_obs = kwargs['readcols'][0], kwargs['readcols'][1]
        wave,flux=np.loadtxt(filename,skiprows=1,usecols=(col_wave_obs, col_flux_obs),unpack=True,dtype=float)
    if 'error' in kwargs:
        col_error = int(kwargs['error'])
        error = np.loadtxt(filename,skiprows=1,usecols=(col_error),unpack=True,dtype=float)
        return wave,flux,error
    else:
        return wave,flux


def main():
    description = '''\
description:

  This scipt can be used for:
   - normalization
   - relative flux calibration
     (this req. a calibrated comparison spectrum;
      use -c comp.spec and -fc res)

  left/right-click to select/deselect continuum points
  hit enter to fit polynomial through continuum points
  hit p to save continuum points to file
  hit i to read continuum points from file
  hit r to reset to the original spectrum
  hit n to normalize (or show instr. profile in fluxcal mode)
  hit w to write to file

  you can use the standard mpl window buttons to zoom in
  if you want to precisely place continuum points

  in case of flux calibration, make sure that:
   - there is no radial velocity shift between spectra
   - and that the comparison resolution matches the observed spectrum
    '''
    epilog = '''\
examples:
  3_fluxcal.py star.dat -c vega_syn_kurucz.txt -fc 7 -w 7
  3_fluxcal.py star.dat -c /home/castor/scripts/evolved/vega_syn_kurucz.txt -fc 7 -w 7
  /home/castor/scripts/evolved/3_fluxcal.py star.dat -c /home/castor/scripts/evolved/vega_syn_kurucz.txt -fc 7 -w 7
  3_fluxcal.py star.dat -c vega_syn_kurucz.txt -fc 7 -w 7 --fix
    '''
  #fluxcal.py HI.20070504.38715_1_01_flux.txt
  #fluxcal.py HI.20070504.38715_1_01_flux.txt -w 0.2
  #fluxcal.py HI.20070504.38715_1_01_flux.txt --fix
  #fluxcal.py HI.20070504.38715_1_01_flux.txt -c HI.20050810.20686_flux.nspec.txt --fix -err 2
  #fluxcal.py mars_20s_0_stacked_sky_1d.dat -c Kurucz_sun2011.txt -fc 3 -w 10 -co 1 2 -scw 10

    parser = argparse.ArgumentParser(
          prog='fluxcal.py',
          formatter_class=argparse.RawDescriptionHelpFormatter,
          description=description,
          epilog=epilog)

    parser.add_argument("--flat",
                        nargs=1, default=False,
                        help="Use flat-fielded flux (specify column number)")
    parser.add_argument("-d","--degree",
                        nargs="?", default=3,type=int,
                        help="Polynomial degree; [1,3,5] (default: 3)")
    parser.add_argument("--fix",
                        action='store_true', default=False,
                        help="Use fixed points instead of continuum median")
    parser.add_argument("-w","--window",
                        nargs="?", default=1,type=float,
                        help="Window size for computing median in Angstrom (default: 1)")
    parser.add_argument("-c","--comparison",
                        nargs="?", default=False,
                        help="Also plot a (normalized) comparison spectrum")
    parser.add_argument("-vrad",
                        nargs=1, default=False,type=float,
                        help="Shift comp. spectrum by vrad")
    parser.add_argument("-err","--error",
                        nargs=1, default=False,
                        help="Also normalize and save error column (specify column number)")
    parser.add_argument("-fc","--fluxcal",
                        nargs="?", default=False,
                        help="Fit spline to comp. spectrum and save continuum/spline_comp, give instr. dlam")
    parser.add_argument("-scw","--scalecompwave",
                        nargs="?", default=False,
                        help="Multiply comparison wave by float")

    parser.add_argument('-co',"--cols_obs",
                    nargs=2, default=False,
                    help="Set col. number for obs. WAVE FLUX (e.g. 1,2)")

    parser.add_argument('argv', nargs='*')
    args = parser.parse_args()
    print("args:\n",args)
    filename = args.argv[0]

    read_kwargs = {}
    if args.flat:
        read_kwargs.update({'flat':args.flat[0]})
    if args.cols_obs:
        read_kwargs.update({'readcols':[int(args.cols_obs[0]),int(args.cols_obs[1])]})
    if args.error:
        read_kwargs.update({'error':args.error[0]})
        wave,flux,error=read_spec(filename,**read_kwargs)
    else:
        wave,flux=read_spec(filename,**read_kwargs)

    kwargs_ontype = {}
    if args.comparison:
        filename_comp = args.comparison
        print("> reading comparison spectrum")
        wave_comp,flux_comp = np.loadtxt(filename_comp,skiprows=1,usecols=(0,1),
                                     unpack=True,dtype=float)
        if args.scalecompwave:
            wave_comp = wave_comp*float(args.scalecompwave)
        if args.fluxcal:
            dlam = float(args.fluxcal)
            if dlam > 0.001:
                midwave_comp = (wave_comp[0] + wave_comp[-1])/2
                print("> comparsion central wavelength",midwave_comp)
                res_corr_comp = midwave_comp/dlam
                print("> convolving comparison with gaussian; this might take a few minutes")
                wave_comp,flux_comp = convolve(wave_comp,flux_comp, res_corr_comp)
        if args.vrad:
            wave_comp = shift_to_obs(wave_comp,args.vrad[0])
            print("> shifted comparison to obs by",args.vrad[0],"km/s")
        kwargs_ontype.update({'comp':[wave_comp,flux_comp]})

    spectrum, = plt.plot(wave,flux,'k-',
                         lw=0.5,marker='.',ms=0.2,label='spectrum')
    plt.title(filename[:20])

    mask = np.where( (np.isfinite(wave) == True) & (np.isfinite(flux) == True) )
    if args.error:
        mask = np.where((np.isfinite(wave) == True) &
                        (np.isfinite(flux) == True) &
                        (np.isfinite(error) == True))
        error = error[mask]
    wave,flux = wave[mask],flux[mask]
    xmin = min(wave[mask])
    xmax = max(wave[mask])
    xran = xmax - xmin
    ymin = min(flux[mask])
    ymax = max(flux[mask])
    if args.flat:
        cut = int(len(flux[mask])*0.993)
        ymax = max(np.sort(flux[mask])[:cut])
    yran = ymax - ymin

    plt.ylim(ymin-0.01*yran,ymax+0.01*yran)

    # Connect the different functions to the different events
    global clicked_points
    clicked_points = []
    plt.gcf().canvas.mpl_connect('key_press_event',lambda event: 
ontype(event,wave,flux,ymin,ymax,yran,filename,vars(args),kwargs_ontype))
    plt.gcf().canvas.mpl_connect('button_press_event',lambda event: onclick(event,wave,flux,**vars(args)))
    plt.gcf().canvas.mpl_connect('pick_event',lambda event: onpick(event,**vars(args)))
    plt.tight_layout()
    plt.show() # show the window

if __name__ == "__main__":
    main()
