#!/usr/bin/env python3

def read_spec(fname,**kwargs):
    import pandas as pd
    " read synspec output spectrum\
      or any spectrum with two columns: wave flux "
    df0  =   pd.read_table(fname,sep='\s+',header=None,comment='#')
    col_wave_obs, col_flux_obs = 0, 1
    if 'readcols' in kwargs:
        col_wave_obs, col_flux_obs = kwargs['readcols'][0], kwargs['readcols'][1]
    rename_dict = {col_wave_obs: 'WAVE', col_flux_obs: 'FLUX'}
    df0 = df0.rename(rename_dict, axis='columns')
    df0[df0['FLUX'] < 0] = float('NaN')
    df0 = df0.dropna()
    return df0

def write_synspec(df_out,savename):
    print('saving as:\n',savename)
    df_out.to_csv(savename,sep=' ', float_format = '%.7e',
                  index=False, header=False)

def est_SNR(flux,**kwargs):
    ''' this is mostly based on https://stdatu.stsci.edu/vodocs/der_snr.pdf \
        works only if there are no gaps in spectrum'''
    # Values that are equal or below zero (padded) are skipped
    import numpy as np
    import scipy.special as sp
    #flux = np.array(flux[np.where(flux > 0.0)])
    n = len(flux)
    # For spectra shorter than this, no value can be returned
    order = 3
    npixmin = 4
    if 'order' in kwargs:
        order = kwargs['order']
    if order == 4: npixmin = 6
    if (n>npixmin):
        signal = np.median(flux)
        if   order == 4:
            # using 4th order of median (tracks spectral lines better, may overestimate SNR)
            # testing with gaussian noise indicates this is more reliable at least in UV and SNR=30
            f4 = 1.482602 / np.sqrt(sp.binom(3, 0)**2 + sp.binom(3, 1)**2 + sp.binom(3, 2)**2 + sp.binom(3, 3)**2)
            noise  = f4 * np.median(np.abs(3.0 * flux[4:n-4] - flux[2:n-6] - 3.0 * flux[6:n-2] + flux[8:n]))
        elif order == 3:
            # using 3rd order of median (linear slopes not counted as noise, may underestimate SNR)
            f3 = 1.482602 / np.sqrt(sp.binom(2, 0)**2 + sp.binom(2, 1)**2 + sp.binom(2, 2)**2)
            noise  = f3 * np.median(np.abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
        return float(signal / noise), noise
    else:
        return float('NaN'), float('NaN')
    
def estimate_continuum(wave, flux, degree):
    from scipy.interpolate import splrep,splev
    import numpy as np
    from astropy.stats import sigma_clip, mad_std

    inilen = len(wave)
    wave, flux = np.array(wave), np.array(flux)
    
    contapprox = convolve_box(flux,int(len(flux)/2))
    
    flux_pn = flux/contapprox
    
    fmedian = np.median(flux_pn)
    mask = np.where((flux_pn > fmedian))
    wave_edit, flux_edit = wave[mask], flux[mask]
    
    
    filtered_data = sigma_clip(flux_edit, sigma=3, stdfunc=mad_std)
    wave_edit, flux_edit = wave_edit[filtered_data.mask], flux_edit[filtered_data.mask]
    
    print("len(wave),len(wave)/inilen\n",len(wave),len(wave)/inilen)

    #flux_syn = convolve_box(flux_syn,npix)
    #error_obs, SNR = est_SNR_array(wave, flux, 1)
    
    x = [it for idx,it in enumerate(wave) if idx % 100 == 0]
    y = [it for idx,it in enumerate(flux) if idx % 100 == 0]
    #x = wave_edit
    #y = flux_edit
    print(x)
    print("len(x),len(y)",len(x),len(y))
    spline = splrep(x,y,k=degree)

    continuum = splev(wave,spline)
    print(len(continuum))
    return continuum, x, y, contapprox

def est_SNR_array(wave_obs, flux_obs, winsizeaa):
    import numpy as np
    """
    Estimate SNR/uncertainty for each flux point by applying
    est_SNR in window around each point
    winsizeaa: window size in Angstrom
    """
    # estimate SNR and errors for each wavelength/flux
    avspacing = np.median([wave_obs[i+1] - wave_obs[i] for i in range(0,len(wave_obs)-1)])
    winsize = max(min(int(winsizeaa/avspacing),len(wave_obs)/10),20)
    flux_obs = np.array(flux_obs)
    SNR_list = [est_SNR(flux_obs[i:i+winsize],**{'order':3}) for i in range(0,len(flux_obs))]
    SNR = np.array([item[0] for item in SNR_list])
    error_obs = np.array([item[1] for item in SNR_list])
    # replace nan values by interpolated values
    nans, x= np.isnan(error_obs), lambda z: z.nonzero()[0]
    error_obs[nans]= np.interp(x(nans), x(~nans), error_obs[~nans])
    return error_obs, SNR

def convolve_FWHM(wavec, fluxc, res):
    ''' Apply instr. convolution with constant FWHM \
        FWHM = ((wave[-1]+wave[0])/2) / res'''
    import numpy as np
    from scipy import interpolate
    from PyAstronomy import pyasl
    # interpolate to equidistant steps using spline
    midwave_syn = (wavec[0] + wavec[-1])/2
    dlam =    midwave_syn / res
    
    if len(wavec)>1:
        stepc   =   abs(wavec[1]-wavec[0])
        if stepc == 0.0:
            stepc =   abs(wavec[2]-wavec[1])
        if dlam < 50:
            stepc = dlam/200
        #print("stepc=",stepc)
        xnew    =   np.arange(wavec[0], wavec[-1], stepc)
        ynew = np.interp(xnew, wavec, fluxc)
        maxsig = 10
        r2, fwhm = pyasl.instrBroadGaussFast(xnew, ynew, res,
          edgeHandling="firstlast", fullout=True, maxsig=maxsig)
        form5   = lambda x: "%.5f" % x
        #print("Gaussian kernel FWHM:", form5(fwhm), 'A')
        step_final = stepc
        if dlam < 50:
            step_final = dlam/5
        #print("step_final=",step_final)
        xnew_app = np.arange(xnew[0], xnew[-1], step_final)
        ynew_app = np.interp(xnew_app, xnew, r2)
        return xnew_app, ynew_app
    else:
        print("Warning: Empty range during convolution")
        return wavec, fluxc

def convolve_box(signal,npix):
    from astropy.convolution import convolve, Box1DKernel
    # preserve_nan=True -> don't set nan values to 0 after conv., instead to nan
    # could also interpolate by setting nan_treatment = 'interpolate'
    smoothed_signal = convolve(signal, Box1DKernel(npix), preserve_nan=True)
    return smoothed_signal

class _Gdl:
  import six.moves as smo
  import numpy as np
  def __init__(self, vsini, epsilon):
    """
      Calculate the broadening profile.
      
      Parameters
      ----------
      vsini : float
          Projected rotation speed of the star [km/s]
      epsilon : float
          Linear limb-darkening coefficient
    """
    self.vc = vsini / 299792.458
    self.eps = epsilon
  
  def gdl(self, dl, refwvl, dwl):
    import numpy as np

    """
      Calculates the broadening profile.
      
      Parameters
      ----------
      dl : array
          'Delta wavelength': The distance to the reference point in
          wavelength space [A].
      refwvl : array
          The reference wavelength [A].
      dwl : float
          The wavelength bin size [A].
      
      Returns
      -------
      Broadening profile : array
          The broadening profile according to Gray. 
    """
    self.dlmax = self.vc * refwvl
    self.c1 = 2.*(1.- self.eps) / (np.pi * self.dlmax * (1. - self.eps/3.))
    self.c2 = self.eps / (2.* self.dlmax * (1. - self.eps/3.))
    result = np.zeros(len(dl))
    x = dl/self.dlmax
    indi = np.where(np.abs(x) < 1.0)[0]
    result[indi] = self.c1*np.sqrt(1. - x[indi]**2) + self.c2*(1. - x[indi]**2)
    # Correct the normalization for numeric accuracy
    # The integral of the function is normalized, however, especially in the case
    # of mild broadening (compared to the wavelength resolution), the discrete
    # broadening profile may no longer be normalized, which leads to a shift of
    # the output spectrum, if not accounted for.
    result /= (np.sum(result) * dwl)
    return result

def rotBroad(wvl, flux, epsilon, vsini, edgeHandling="firstlast"):
  import six.moves as smo
  import numpy as np
  """Apply rotational broadening to a spectrum."""
  # Check whether wavelength array is evenly spaced
  sp = wvl[1::] - wvl[0:-1]
  if abs(max(sp) - min(sp)) > 1e-6:
    print("Input wavelength array is not evenly spaced")
  if vsini <= 0.0:
    print("vsini must be positive.")
  if (epsilon < 0) or (epsilon > 1.0):
    print("Linear limb-darkening coefficient, epsilon, should be '0 < epsilon < 1'")
  
  # Wavelength binsize
  dwl = wvl[1] - wvl[0]
  
  # Indices of the flux array to be returned
  validIndices = None
  
  if edgeHandling == "firstlast":
    # Number of bins additionally needed at the edges 
    binnu = int(np.floor(((vsini / 299792.458) * max(wvl)) / dwl)) + 1
    # Defined 'valid' indices to be returned
    validIndices = np.arange(len(flux)) + binnu
    # Adapt flux array
    front = np.ones(binnu) * flux[0]
    end = np.ones(binnu) * flux[-1]
    flux = np.concatenate( (front, flux, end) )
    # Adapt wavelength array
    front = (wvl[0] - (np.arange(binnu) + 1) * dwl)[::-1]
    end = wvl[-1] + (np.arange(binnu) + 1) * dwl
    wvl = np.concatenate( (front, wvl, end) )
  elif edgeHandling == "None":
    validIndices = np.arange(len(flux))
  else:
    print("Edge handling method not supported")
  
  result = np.zeros(len(flux))
  gdl = _Gdl(vsini, epsilon)
  
  considerpix = int((vsini/10)*9)
  if considerpix < 2: considerpix = 2
  for i in smo.range(considerpix+1):
    dl = (wvl[i] - wvl[:i+considerpix])
    #print('dl:',dl)
    #print('wvl[i]:',wvl[i])
    #print('dwl:',dwl)
    g = gdl.gdl(dl, wvl[i], dwl)
    result[i] = np.sum(flux[:i+considerpix] * g)
  for i in smo.range(considerpix+1,len(flux)-(considerpix+1)):
    dl = (wvl[i] - wvl[(i-considerpix):(i+considerpix)])
    if i%5000 == 0: print('wvl[i]:',wvl[i])
    g = gdl.gdl(dl, wvl[i], dwl)
    result[i] = np.sum(flux[(i-considerpix):(i+considerpix)] * g)
  for i in smo.range(len(flux)-(considerpix+1),len(flux)):
    dl = (wvl[i] - wvl[(i-considerpix):])
    g = gdl.gdl(dl, wvl[i], dwl)
    result[i] = np.sum(flux[(i-considerpix):] * g)
  result *= dwl
  return result[validIndices]

def apply_rotbroad_slow(wavec, fluxc, vsini, epsilon, wave_obs):
    '''accounts for doppler shift wave dependency, but is super slow'''
    from PyAstronomy import pyasl
    from scipy import interpolate
    import numpy as np
    print('applying wavelength-dependent rotational broadening')
    print('vsini =',vsini)
    stepc   =   abs(wavec[1]-wavec[0])
    if stepc == 0.0:
        stepc =   abs(wavec[2]-wavec[1])
    lwobs = int(len(wave_obs)*(1/10))
    # Save time by not oversampling too much
    stepc = abs(wave_obs[lwobs+1]-wave_obs[lwobs])/2
    print("Wavelength step:",stepc)
    tck     =   interpolate.splrep(wavec, fluxc, s=0)
    xnew    =   np.arange(wavec[0], wavec[-1], stepc)
    ynew    =   interpolate.splev(xnew, tck, der=0)
    flux_conv = rotBroad(xnew, ynew, epsilon, vsini, edgeHandling='firstlast')
    return xnew, flux_conv

def apply_rotbroad_fast(wavec, fluxc, vsini, epsilon, effWvl):
    '''neglects doppler shift wave dependency, but is super fast'''
    from PyAstronomy import pyasl
    from scipy import interpolate
    import numpy as np
    print('applying wave-indep. rot. broadening')
    print('vsini =',vsini)
    print('effWvl =',effWvl)
    stepc   =   abs(wavec[1]-wavec[0])
    if stepc == 0.0:
        stepc =   abs(wavec[2]-wavec[1])
    tck     =   interpolate.splrep(wavec, fluxc, s=0)
    xnew    =   np.arange(wavec[0], wavec[-1], stepc)
    ynew    =   interpolate.splev(xnew, tck, der=0)
    flux_conv = pyasl.fastRotBroad(xnew, ynew, epsilon, vsini, effWvl=effWvl)
    return xnew, flux_conv

def edit_sec(wave_edit,flux_edit,nsteps,anotherfunc, extraArgs):
    ''' apply different convolution for single panels \
        need to improve the ranges in which the FWHM is constant '''
    import numpy as np
    wave_start0 = wave_edit[0]
    wave_end0   = wave_edit[-1]
    wave_range  = wave_end0 - wave_start0
    wave_step0  = wave_range/nsteps
    fullwave = []
    fullflux = []
    for i in range(0,nsteps):
        wave_start  = wave_start0 + i      * wave_step0
        wave_end    = wave_start0 + (i+1)  * wave_step0
        cut_to_panel = np.where((wave_start<wave_edit)&(wave_edit<wave_end))
        print('Conv. range:',round(wave_start,2),'-',round(wave_end,2))
        wave_edit_cut = wave_edit[cut_to_panel]
        flux_edit_cut = flux_edit[cut_to_panel]

        wave_conv,flux_conv = anotherfunc(wave_edit_cut,
                                          flux_edit_cut,
                                          *extraArgs)

        cut_conv = np.where((wave_start<wave_conv)&(wave_conv<wave_end))
        fullwave.append(wave_conv[cut_conv])
        fullflux.append(flux_conv[cut_conv])
    fullwave = [item for sublist in fullwave for item in sublist]
    fullflux = [item for sublist in fullflux for item in sublist]
    return np.array(fullwave), np.array(fullflux)

def shift_to_obs(waves,vrad):
    import numpy as np
    c = 299792.458 # km/s
    fac_toobs = np.sqrt((1+(vrad/c))/(1-(vrad/c)))
    shifted = np.multiply(waves,fac_toobs)
    return shifted


def main():
    import pandas as pd
    import numpy as np
    import argparse,textwrap,time,os
    time_start = time.time()
    parser = argparse.ArgumentParser(
         prog='edit_syn7.py',
         formatter_class=argparse.RawDescriptionHelpFormatter,
         epilog=textwrap.dedent('''\
            description:
                 edit input spectrum: e.g.
                   apply rot. broadening
                   apply instr. broadening
                   estimate SNR
                   mulitply/divide by second spectrum
                   mulitply flux by constant
            examples:
                 edit_syn7.py -spec syn.7 -res 6
         '''))
    parser.add_argument('-test',
                        action='store_true', default=False,
                        help="Testing purpose")
    parser.add_argument('-esterr',
                        action='store_true', default=False,
                        help="(Over)write error column with error estimated from flux and wave arrays")
    parser.add_argument('-c','--comparison',
                        nargs=1, default=False,
                        help="Comparsion spectrum, e.g. obs")
    parser.add_argument('-spec',
                        nargs=1, default=False,
                        help="Input spectrum")
    parser.add_argument('-vsini',
                        nargs='?', default=False,
                        help="Apply rotational broadening")
    parser.add_argument('-res',
                        nargs='?', default=False,
                        help="Apply instrumental broadening")
    parser.add_argument('-res12',
                        nargs=2, default=False,
                        help="Apply instrumental broadening to already broadened spectrum\
                              supply target and current resolution R")

    parser.add_argument('-sb',"--smooth_box",
                        nargs=1, default=[0],type=float,
                        help="Smooth spectra using box filter of specified width in pixels")
    parser.add_argument('-reverse',
                        action='store_true', default=False,
                        help="Reverse order of wave and flux")
    parser.add_argument('-wc','--wave_const',
                        nargs=1, default=False,
                        help="Multipy input wave by constant")
    parser.add_argument('-fc','--flux_const',
                        nargs=1, default=False,
                        help="Multipy input flux by constant")
    parser.add_argument('-m','--multiply',
                        nargs=1, default=False,
                        help="Multipy input spectrum with other spectrum, e.g. ISM extinction")
    parser.add_argument('-d','--divide',
                        nargs=1, default=False,
                        help="Divide input spectrum by other spectrum, e.g. instrumental profile")
    parser.add_argument('-toHlam',
                        action='store_true', default=False,
                        help="Convert flux units from ? to Hlam? (mult. by 4pi)")
    parser.add_argument('-reg',
                        action='store_true', default=False,
                        help="interpolate to regular grid (using av. distance)")
    parser.add_argument('-aa_to_ev',
                        action='store_true', default=False,
                        help="Convert wavelength in Angstrom to energy in eV")
    parser.add_argument('-cols1',
                    nargs=2, default=False,
                    help="Set col. numbers for input spectrum WAVE FLUX (e.g. 1 2)")
    parser.add_argument('-cols2',
                    nargs=2, default=False,
                    help="Set col. numbers for second spectrum WAVE FLUX (e.g. 1 2)")
    parser.add_argument('argv', nargs='*')
    args = parser.parse_args()
    print('args\n:',args)

    specname = args.spec[0]
    read_kwargs = {}
    if args.cols1:
        read_kwargs.update({'readcols':[int(args.cols1[0]),int(args.cols1[1])]})
    spectrum = read_spec(specname,**read_kwargs)
    print(spectrum)
    wave_syn = np.array(spectrum['WAVE'])
    flux_syn = np.array(spectrum['FLUX'])

    sname_ext = ''

    if args.comparison:
        specname2 = args.comparison[0]
        wave_obs, flux_obs, err_obs  = np.loadtxt(specname2,
                               usecols=(0,1,2),
                               unpack=True,
                               dtype=float)
        wave_obs, flux_obs, err_obs = np.array(wave_obs), np.array(flux_obs), np.array(err_obs)

    if args.test:
        flux_syn = np.interp(wave_obs, wave_syn, flux_syn)
        wave_syn, flux_syn = np.array(wave_syn), np.array(flux_syn)
        matched_wave,matched_data,matched_error,matched_model = \
        match_data_models( wave_obs, flux_obs, np.ones(len(wave_obs)), err_obs, wave_obs, flux_syn, max(min(wave_syn),min(wave_obs)),min(max(wave_syn),max(wave_obs)),
                         saveDowngradedModel = True, downgradedModelFile = "DGmodel.txt")

    if args.wave_const:
        wave_syn  = np.multiply(wave_syn,float(args.wave_const[0]))
        sname_ext = sname_ext + '.wavec'
    if args.flux_const:
        flux_syn  = flux_syn * float(args.flux_const[0])
        sname_ext = sname_ext + '.fluxc'
    if args.toHlam:
        # actually F_lam = H_lam *4pi, so this converts H_lam to F_lam
        flux_syn = flux_syn * 4 * 3.14159265359
        sname_ext = sname_ext + '.Hlam'
    if args.aa_to_ev:
        v_light = 299792458 #m/s
        wave_syn = v_light/(wave_syn*1e-10)*4.135667662e-15
        sname_ext = sname_ext + '.ev'

    if args.res:
        sname_ext = sname_ext + '.R' + str(args.res)
        resolution = float(args.res)
        if resolution > 50:
            resolution = int(resolution)
            wave_syn,flux_syn = edit_sec(wave_syn,flux_syn,100,convolve_FWHM,[resolution])
        else:
            dlam = resolution
            midwave_syn = (wave_syn[0] + wave_syn[-1])/2
            res_corr_syn = midwave_syn/dlam
            wave_syn,flux_syn = convolve_FWHM(wave_syn,flux_syn, res_corr_syn)

    if args.res12:
        ress = sorted([float(i) for i in args.res12])
        res_tar = ress[0]
        res_ini = ress[1]
        res_apply = 1/(1/res_tar - 1/res_ini)
        print("Convolving with R=",res_apply)
        wave_syn,flux_syn = edit_sec(wave_syn,flux_syn,100,convolve_FWHM,[res_apply])
        sname_ext = sname_ext + '.R' + str(res_tar)

    if args.smooth_box[0] >= 2:
        npix = int(args.smooth_box[0])
        flux_syn = convolve_box(flux_syn,npix)
        sname_ext = sname_ext + '.sb'+str(args.smooth_box[0])

    if args.reverse:
        wave_syn = np.flip(wave_syn)
        flux_syn = np.flip(flux_syn)
        sname_ext = sname_ext + '.reverse'
    if args.multiply:
        specname2 = args.multiply[0]
        print("Mulitplying with provided spectrum",specname2)
        read_kwargs = {}
        if args.cols2:
            read_kwargs.update({'readcols':[int(args.cols2[0]),int(args.cols2[1])]})
        spectrum2 = read_spec(specname2,**read_kwargs)
        spectrum2 = spectrum2.sort_values(by=['WAVE'])
        wave_syn2 =  np.array(spectrum2['WAVE'])
        flux_syn2  = np.array(spectrum2['FLUX'])
        flux_syn2_intp = np.interp(wave_syn, wave_syn2, flux_syn2)
        flux_syn = flux_syn*flux_syn2_intp
        sname_ext = sname_ext + '.multiply'
    if args.divide:
        specname2 = args.divide[0]
        print("Mulitplying with provided spectrum",specname2)
        if args.cols2:
            read_kwargs.update({'readcols':[int(args.cols2[0]),int(args.cols2[1])]})
        spectrum2 = read_spec(specname2,**read_kwargs)
        spectrum2 = spectrum2.sort_values(by=['WAVE'])
        wave_syn2 =  np.array(spectrum2['WAVE'])
        flux_syn2  = np.array(spectrum2['FLUX'])
        flux_syn2_intp = np.interp(wave_syn, wave_syn2, flux_syn2)
        flux_syn = flux_syn/flux_syn2_intp
        sname_ext = sname_ext + '.divide'
    if args.reg:
        wave_syn_reg = np.linspace(min(wave_syn), max(wave_syn), int(len(wave_syn)*1.1))
        flux_syn_reg = np.interp(wave_syn_reg, wave_syn, flux_syn)
        wave_syn = wave_syn_reg
        flux_syn = flux_syn_reg
        sname_ext = sname_ext + '.reg'

    if args.esterr:
        sname_ext = sname_ext + '.esterr'
        error_syn, SNR = est_SNR_array(wave_syn, flux_syn, 1)

    savename = os.path.splitext(specname)[0] + sname_ext
    #d = {'WAVE': wave_syn, 'FLUX': flux_syn}
    np.savetxt(savename,np.array([wave_syn, flux_syn]).T)
    #if args.esterr:
    #    d = {'WAVE': wave_syn, 'FLUX': flux_syn, 'ERROR': error_syn}
    #outframe = pd.DataFrame(data=d)
    #write_synspec(outframe,savename)

if __name__ == "__main__":
    main()

