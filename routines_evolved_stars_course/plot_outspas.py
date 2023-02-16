#!/usr/bin/env python3

def norm_lin(wave,flux,wl1,wl2,wl3,wl4):
    import numpy as np
    rlow    = np.where((wl1 <= wave) & (wave <= wl2))
    rhigh   = np.where((wl3 <= wave) & (wave <= wl4))
    flux_s1 = flux[rlow]
    flux_s2 = flux[rhigh]
    
    if wl1 < 3000:
        cut_max = 0.9
    if wl1 >= 3000:
        cut_max = 1/3
    
    flux_s1_sc = np.sort(flux_s1)[int(len(flux_s1)*cut_max):]
    flux_s2_sc = np.sort(flux_s2)[int(len(flux_s2)*cut_max):]
    
    y  = ((np.average(flux_s2_sc) - np.average(flux_s1_sc)) /\
        (wl3/2 + wl4/2 - wl1)) \
        * (wave - wl1) + np.average(flux_s1_sc)
    return y

def convolve(wavec, fluxc, res):
    from scipy import interpolate
    from PyAstronomy import pyasl
    form5   = lambda x: "%.5f" % x
    # interpolate to equidistant steps using spline
    stepc   =   abs(wavec[1]-wavec[0])
    tck     =   interpolate.splrep(wavec, fluxc, s=0)
    xnew    =   np.arange(wavec[0], wavec[-1], stepc)
    ynew    =   interpolate.splev(xnew, tck, der=0)
    # apply instr. Gauss
    maxsig = 50
    r2, fwhm = pyasl.instrBroadGaussFast(xnew, ynew, res,
          edgeHandling="firstlast", fullout=True, maxsig=maxsig)
    print("Gaussian kernel FWHM:", form5(fwhm), 'A')
    return xnew, r2

# apply different convolution for single panels
# need to improve the ranges in which the FWHM is constant
def edit_sec(wave_edit,flux_edit,resolution):
    form2   = lambda x: "%.2f" % x
    #cut_to_obs  = np.where((cut_wl-0.01<=wave_edit)&(wave_edit+0.01<=end_wl))
    #wave_edit   = wave_edit[cut_to_obs]
    wave_start0 = wave_edit[0]
    wave_end0   = wave_edit[-1]
    wave_range  = wave_end0 - wave_start0
    wave_step0  = wave_range/100
    fullwave = []
    fullflux = []
    for i in range(0,100):
        wave_start  = wave_start0 + i      * wave_step0
        wave_end    = wave_start0 + (i+1)  * wave_step0
        cut_to_panel = np.where((wave_start<wave_edit)&(wave_edit<wave_end))
        print('Conv. range:',form2(wave_start),'-',form2(wave_end))
        wave_conv,flux_conv = \
            convolve(wave_edit[cut_to_panel],
                     flux_edit[cut_to_panel],
                     resolution)
        cut_conv = np.where((wave_start<wave_conv)&(wave_conv<wave_end))
        fullwave.append(wave_conv[cut_conv])
        fullflux.append(flux_conv[cut_conv])
    fullwave = [item for sublist in fullwave for item in sublist]
    fullflux = [item for sublist in fullflux for item in sublist]
    return np.array(fullwave), np.array(fullflux)

def multiple_replace(string, rep_dict):
    import re
    pattern = \
        re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

def read_rpar(nrows,spasout):
    import pandas as pd
    rindex, x0index,parlist  =  [], [], []
    lookup0, lookup1 = 'x0=', '# wave'
    elements = []
    wave_cen = []
    Rs       = []
    with open(spasout) as myFile:
        for index, line in enumerate(myFile, 1):
            if lookup1 in line:
                rindex.append(index)
            if lookup0 in line:
                x0index.append(index)
                edit0 = line.split('\t')
                rdict = {'x0=':'',
                         'R=':'',
                         '\n':'',
                         '# ':'',
                         'Continuum=':''}
                edit = [multiple_replace(item,rdict) for item in edit0]
                edit = [item.split(',') for item in edit]
                edit = [item for sublist in edit for item in sublist]
                elements.append(edit[0])
                wave_cen.append(float(edit[1]))
                Rs.append(float(edit[2]))
                parlist.append(edit)
    c = [list(item) for item in list(zip(x0index,rindex))]
    index_start = []
    index_stop  = []
    for i in range(0,len(c)):
        if i < len(c)-1:
            index_start.append(c[i][1])
            index_stop.append(c[i+1][0]-1)
        else:
            index_start.append(c[i][1])
            index_stop.append(nrows)
    ranges = [list(item) for item in list(zip(index_start,index_stop))]
    d = [list(item) for item in list(zip(ranges,parlist))]
    
    dataframe = pd.DataFrame({'index_start': index_start,
                       'index_stop': index_stop,
                       'element': elements,
                       'wave_cen': wave_cen,
                       'R': Rs})
    
    return d, dataframe

def rawincount(filename):
    from itertools import (takewhile,repeat)
    f = open(filename, 'rb')
    bufgen = \
        takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

def read_ranges(data,par):
    import numpy as np
    waves,flux_obss,flux_models = [],[],[]
    for i in range(0,len(par)):
        wave,flux_obs,flux_model  = np.loadtxt(data[par[i][0][0]:par[i][0][1]],
                                               usecols = (0,1,2),
                                               unpack=True)
        waves.append(wave)
        flux_obss.append(flux_obs)
        flux_models.append(flux_model)
    return waves,flux_obss,flux_models

def meta_data(Teff,logg,logy,vrad,vrot,zeta,data_edit):
    import matplotlib.pyplot as plt
    ''' create and plot a LaTex table with metadata '''
    import matplotlib as mpl
    #mpl.rc('text', usetex=True)
    form1   = lambda x: "%.1f" % x
    fig = plt.figure()
    star = 'Teststar'
    obsname = 'Testobs'
    basename_0 = 'Testmodel'
    resolution = 'Testres'
    table = r'{\renewcommand{\arraystretch}{1.3}'\
    + r'\begin{table} \begin{tabular}{|l|l|} \hline ' \
    + r'Star & '+multiple_replace(star,{'_':r'\underline{{ }}','-':'$-$'})+r' \\ \hline ' \
    + r'Observation & '+multiple_replace(obsname,{'_':r'\underline{{ }}','-':'$-$'})+r' \\ \hline ' \
    + r'Model & '+multiple_replace(basename_0,{'_':r'\underline{{ }}','-':'$-$'})+r' \\ \hline ' \
    + r'$v$$_\mathrm{rad}$ (km/s) & '+str(form1(vrad))+r'  \\ \hline ' \
    + r'$R$ & '+str(resolution)+r'  \\ \hline '
    table = table + r'\end{tabular} \end{table}}'
    print(table)
    data_string_1 = ''.join(data_edit[:3])
    data_string_2 = ''.join(data_edit[3:])
    table = data_string_1 + '\n' + data_string_2
    ax = fig.add_subplot(111)
    ax.text(0.5, 0.5,table,
            ha='center',va='center',size=12)
    ax.axis('off')
    plt.tight_layout(pad=0)
    plt.close()
    #mpl.rc('text', usetex=False)
    return fig

def plot_singlecolumn(dataframe,**kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    wave = kwargs.get('wave', [])
    flux = kwargs.get('flux', [])

    replace_dict = {'Ha':r'H$\alpha$/HeII',
                    'Hb':r'H$\beta$/HeII',
                    'Hg':r'H$\gamma$/HeII',
                    'Hd':r'H$\delta$/HeII',
                    'He':r'H$\epsilon$/HeII',
                    'Ha + HeI':r'H$\alpha$/HeII + HeI',
                    'Hb + HeI':r'H$\beta$/HeII + HeI',
                    'Hg + HeI':r'H$\gamma$/HeII + HeI',
                    'Hd + HeI':r'H$\delta$/HeII + HeI',
                    'He + HeI':r'H$\epsilon$/HeII + HeI'}
    
    dimx = len(dataframe.values)
    if 'axes' in kwargs:
        ax1 = kwargs['axes']
    else:
        fig =   plt.figure(figsize=(8.,dimx*1.3))
        ax1 =   plt.gca()
        ax1.set_xlabel(r"$\lambda -\lambda _0$ ($\AA$)",
                   color='black',fontsize=13)
        ax1.set_ylabel(r"Normalized Flux + Offset",
                   color='black',fontsize=13)
    alignwave = [dataframe.iloc[k].wave_cen for k in list(range(dimx))]
    if 'center' in kwargs:
        alignwave = [(dataframe.iloc[k].wave[0] + dataframe.iloc[k].wave[-1])/2 for k in list(range(dimx))]
    waveminmin = min([min(dataframe.iloc[k].wave - alignwave[k]) for k in list(range(dimx))])
    wavemaxmax = max([max(dataframe.iloc[k].wave - alignwave[k]) for k in list(range(dimx))])
    xdist = abs(wavemaxmax) + abs(waveminmin)
    #lamrange = 1.1*max(abs(waveminmin),abs(wavemaxmax))
    xmin_panel = -abs(waveminmin)-0.03*xdist
    xmax_panel = abs(wavemaxmax)+0.03*xdist
    ax1.set_xlim(xmin_panel,xmax_panel)

    if 'save_space' in kwargs:
        ypos_tight = []
        for i in range(0,dimx):
            item = dataframe.iloc[i]
            alignwave_ss = item.wave_cen
            if 'center' in kwargs:
                alignwave_ss = (item.wave[0]+item.wave[-1])/2
            xrange_current = item.wave - alignwave_ss
            # "normalizing" by max(item.obs) is not a good idea:
            # e.g. cosmics will lead to a very large max(item.obs)
            # -> very small values for min(item.obs[i]/max(item.obs)
            # still need some way to normalize, since SPAS sometimes has continuum at 2...
            # should only normalize by constant, otherwise wrong continuum placement is corrected ...
            minarray = np.array([min(item.obs[i]/max(item.obs),item.model[i]/max(item.model)) for i in range(0,len(item.obs))])
            ymin = min(min(item.obs/max(item.obs)),min(item.model/max(item.model)))
            ymax = max(max(item.obs/max(item.obs)),max(item.model/max(item.model)))
            #print('ymin,ymax:',ymin,ymax)
            yspace = ymax-ymin
            ybuffer = 0.1
            if i > 0:
                print('\n')
                print('item.wave_cen =',alignwave_ss)
                yspace_prev = ypos_tight[-1]
                item_prev = dataframe.iloc[i-1]
                alignwave_ss_prev = item_prev.wave_cen
                if 'center' in kwargs:
                    alignwave_ss_prev = (item_prev.wave[0]+item_prev.wave[-1])/2
                maxarray_prev = np.array(\
                [max(item_prev.obs[i]/max(item_prev.obs),item_prev.model[i]/max(item_prev.model)) \
                for i in range(0,len(item_prev.obs))])
                xrange_prev = item_prev.wave - alignwave_ss_prev
                try:
                    mask = np.where((min(xrange_current) < xrange_prev) & (xrange_prev < max(xrange_current)))
                    maxarray_prev_int = np.interp(xrange_current,xrange_prev[mask],maxarray_prev[mask])
                except ValueError:
                    maxarray_prev_int = np.interp(xrange_current,xrange_prev,maxarray_prev)
                ydiff_to_prev = maxarray_prev_int - minarray
                ydiff_min = min(ydiff_to_prev)
                ydiff_max = max(ydiff_to_prev)
                #print('ydiff_min =',ydiff_min)
                #print('ydiff_max =',ydiff_max)
                if ydiff_max > 0:
                    yadd =  ydiff_max + ybuffer
                    print('yadd =',yadd)
                    ypos = ypos_tight[-1] + yadd
                else:
                    ypos = ypos_tight[-1] + ybuffer
            else:
                ypos = yspace + ybuffer
            ypos_tight.append(ypos)
        ax1.set_ylim(0,ypos_tight[-1]+ybuffer)
        if 'yrange' in kwargs:
            ax1.set_ylim(0,kwargs['yrange'])
    else:
        offset = max(0.8,1.3/dimx)
        ax1.set_ylim(0,dimx+offset)

    for i in range(0,dimx):
        item = dataframe.iloc[i]
        
        alignwave_p = item.wave_cen
        if 'center' in kwargs:
            alignwave_p = (item.wave[0]+item.wave[-1])/2

        label = item.element
        if label in replace_dict:
            label_new = label.replace(label, replace_dict[label])
        else: label_new = label
        
        if 'save_space' in kwargs:
            txtoffset = 1/20
            txtoffset_0 = 1
            yoffset = ypos_tight[i]
        else:
            txtoffset_0 = 1
            txtoffset = 1/4
            yoffset = i + 1
        
        if not 'axes' in kwargs:
            ax1.plot(item.wave - alignwave_p,
                        item.obs/max(item.model) + yoffset - max(dataframe.iloc[0].model),
                        lw=0.5,color='black')
            ax1.plot(item.wave - alignwave_p,
                        item.model/max(item.model) + yoffset - max(dataframe.iloc[0].model),
                        color='red',lw=0.75,
                        label=label_new)
        else:
            if 'yrange' in kwargs:
                if 'save_space' in kwargs:
                    yshift_std_savespace = kwargs['yrange'] - (ypos_tight[-1] + ybuffer)
                else:
                    yshift_std_savespace = kwargs['yrange'] - (dimx + 1 + offset)
                yoffset = yoffset + yshift_std_savespace
            ax1.plot(item.wave - alignwave_p,
                        item.obs/max(item.model) + yoffset - 1,
                        lw=0.5,color='black')
            ax1.plot(item.wave - alignwave_p,
                        item.model/max(item.model) + yoffset - 1,
                        color='red',lw=0.75,
                        label=label_new)
            if len(wave)>0 and len(flux)>0:
                xmin = min(item.wave)
                xmax = max(item.wave)
                xran = xmax - xmin
                print(xmin,xmax)
                mask = np.where((xmin<wave)&(wave<xmax))
                fnorm = norm_lin(wave,flux,xmin,xmin+0.01*xran,xmax-0.01*xran,xmax)
                fnorm_tospas = norm_lin(item.wave,item.model/max(item.model),
                                        xmin,xmin+0.01*xran,
                                        xmax-0.01*xran,xmax)
                fnorm_tospas = np.interp(wave, item.wave, fnorm_tospas)
                wave_cut = wave[mask]
                flux_norm = flux/fnorm*fnorm_tospas
                flux_plot  = flux_norm[mask]
                ax1.plot(wave_cut - alignwave_p,
                        flux_plot + yoffset - 1,
                        lw=0.5,color='blue')
        plt.annotate(s=label_new,size=13,xy=(0, 1),
                    xycoords='axes fraction',
                    xytext=(xmin_panel+xdist*0.015, yoffset+txtoffset),
                    textcoords='data',
                    ha='left',va='top',color='grey')
        plt.annotate(s=r'$\lambda _0$'+str(round(alignwave_p,1)),
                    size=13,xy=(0, 1),
                    xycoords='axes fraction',
                    xytext=(xmax_panel-xdist*0.015, yoffset+txtoffset),
                    textcoords='data',
                    ha='right',va='top',color='grey')

    ax1_ymin, ax1_ymax = ax1.get_ylim()
    ax1_xmin, ax1_xmax = ax1.get_xlim()
    ax1_xextend = ax1_xmax - ax1_xmin
    ax1_yextend = ax1_ymax - ax1_ymin
    if 'ruler' in kwargs:
        #print('ax1_xmin, ax1_xmax',ax1_xmin, ax1_xmax)
        #print('ax1_ymin, ax1_ymax',ax1_ymin, ax1_ymax)
        y1, y2 = plt.gca().get_window_extent().get_points()[:, 1]
        x1, x2 = plt.gca().get_window_extent().get_points()[:, 0]
        yscale = (y2-y1)/(ax1_ymax-ax1_ymin)
        xscale = (x2-x1)/(ax1_xmax-ax1_xmin)
        #print('xscale,yscale',xscale,yscale)
        fontsize = 12
        # x-extension of one char in x-axis units?!
        xfontsize = fontsize / xscale
        #print('xfontsize',xfontsize)
        mindist = xfontsize * 1.15

        settings = {'color': (0,0,0),
                   'lw': 1,
                   'head_width': 0,
                   'alpha':0.7}
        ruler_cx = 0.7*ax1_xmax
        ruler_cx = ax1_xmin + 0.15*(ax1_xmax-ax1_xmin)
        ruler_xx = xfontsize*2
        ruler_yextend = 1
        while ruler_yextend > ax1_yextend:
            ruler_yextend = ruler_yextend/2
            print(ruler_yextend)
        if 'ruler_length' in kwargs:
            ruler_yextend = kwargs['ruler_length']
        ruler_ystart  = max(int(ax1_ymax)-ruler_yextend,ax1_ymax-np.ceil(ruler_yextend*10)/10-0.1)
        ruler_ystart  = ax1_ymin+0.1
        #ruler_ystart  = ax1_ymax - ruler_yextend - ax1_yextend*0.05
        form2   = lambda x: "%.2f" % x
        form1   = lambda x: "%.1f" % x
        rulerlabel = form2(round(ruler_yextend,2))
        if rulerlabel[-1] == '0':
            rulerlabel =  form1(round(ruler_yextend,1))
        ax1.arrow(ruler_cx, ruler_ystart, 0, ruler_yextend,**settings)
        ax1.arrow(ruler_cx-ruler_xx/2, ruler_ystart, ruler_xx, 0,**settings)
        ax1.arrow(ruler_cx-ruler_xx/2, ruler_ystart+ruler_yextend, ruler_xx, 0,**settings)
        ax1.text(ruler_cx-xfontsize,ruler_ystart+ruler_yextend/2,rulerlabel,
                size=fontsize,rotation=90., ha="center", va="center",alpha=0.7)

    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax1.tick_params(axis='both',which='both',direction='in',
                    top=True,bottom=True, left=True, right=True)
    ax1.tick_params(axis='x',labelsize=12)

    if 'ruler' in kwargs:
        ax1.tick_params(axis='both',which='both',direction='in',
                        top=True,bottom=True, left=True, right=True,labelleft=False)
    if 'nogrid' not in kwargs:
         ax1.xaxis.grid(True,ls='dashed')
    if not 'axes' in kwargs:
        if 'nogrid' not in kwargs:
            ax1.yaxis.grid(True,ls='dashed')
    #ax1.yaxis.set_visible(False)
    #add_scalebar(ax1, matchx=True, matchy=True, hidex=False, hidey=False, **{})
    if not 'axes' in kwargs:
        plt.tight_layout()
        plt.close()
        return fig
    else:
        return ax1_yextend

def plot_multicolumn(dataframe_full,ranges_per_panel,**plot_df_keys):
    ''' create multiple plot_singlecolumn() side by side '''
    import numpy as np
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    dataframe_full = dataframe_full.iloc[::-1]
    nranges = len(dataframe_full.values)
    #ranges_per_panel = 10
    ranges_per_panel = int(nranges/2 + 0.5)
    npanels = int(np.ceil(nranges / ranges_per_panel))
    if 'sortatm' in plot_df_keys:
        # select ranges where element string contains 'HeI' but not 'HeII'
        print("dataframe_full.element:\n",dataframe_full.element)
        mask_HeI = ((dataframe_full.element.str.contains('HeI', regex=True)) &\
                   ~(dataframe_full.element.str.contains('HeII', regex=True)) &\
                   ~(dataframe_full.element.str.contains('He +', regex=True)))
        print("mask_HeI:\n",mask_HeI)
        df_HeI = dataframe_full[mask_HeI]
        df_noHeI = dataframe_full[~mask_HeI]
        dfs_selected = [df_noHeI,df_HeI]
        for item in dfs_selected:
            print(len(item))
        dfs_selected.sort(key=len)
        dfs_selected = dfs_selected[::-1]
        npanels = len(dfs_selected)
    figsize_x = npanels*8.
    figsize_y = ranges_per_panel*1.3
    pname = plt.figure(figsize=(figsize_x,figsize_y))
    outer = gridspec.GridSpec(1, npanels,  wspace=0, hspace=0)
    if 'sortatm' in plot_df_keys:
        for i in range(0,npanels):
            axes_panel = pname.add_subplot(1,npanels,i+1)
            axes_panel.clear()
            plot_df_keys['axes'] = axes_panel
            # only show ruler for first panel
            if i > 0:
                plot_df_keys.pop('ruler', None)
            dataframe_cut = dfs_selected[i]
            yrange = plot_singlecolumn(dataframe_cut,**plot_df_keys)
            plot_df_keys['yrange'] = yrange
            pname.add_subplot(axes_panel)
    else:
        for i in range(0,npanels):
            axes_panel = pname.add_subplot(1,npanels,i+1)
            axes_panel.clear()
            plot_df_keys['axes'] = axes_panel
            if i > 0:
                plot_df_keys.pop('ruler', None)
            dataframe_cut = dataframe_full[i*ranges_per_panel:(i+1)*ranges_per_panel][::-1]
            yrange = plot_singlecolumn(dataframe_cut,**plot_df_keys)
            plot_df_keys['yrange'] = yrange
            pname.add_subplot(axes_panel)
    pname.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='white',which='both',axis='both',color='white')
    #plt.grid(False)
    plt.xlabel("$\lambda -\lambda _0$ ($\AA$)",fontsize=13, labelpad=0)
    plt.ylabel(r"Normalized Flux",fontsize=13, labelpad=0)
    plt.tight_layout(pad=0,w_pad=0.4,h_pad=0)
    #plt.tight_layout()
    plt.close()
    return pname


def plot_singlemosaic(dataframe,ncols,**kwargs):
    import matplotlib.pyplot as plt
    import pandas as pd
    import math
    single = kwargs.get('single', False)
    wave = kwargs.get('wave', [])
    flux = kwargs.get('flux', [])

    replace_dict = {'Ha':r'H$\alpha$/HeII',
                    'Hb':r'H$\beta$/HeII',
                    'Hg':r'H$\gamma$/HeII',
                    'Hd':r'H$\delta$/HeII',
                    'He':r'H$\epsilon$/HeII',
                    'Ha + HeI':r'H$\alpha$/HeII + HeI',
                    'Hb + HeI':r'H$\beta$/HeII + HeI',
                    'Hg + HeI':r'H$\gamma$/HeII + HeI',
                    'Hd + HeI':r'H$\delta$/HeII + HeI',
                    'He + HeI':r'H$\epsilon$/HeII + HeI'}
    
    dimx = len(dataframe.values)
    ncols = 3
    nrows = math.ceil((dimx-ncols)/ncols)+1
    figx = ncols*5.4#6.5#4.8
    figy = nrows*3.6#3.1#3.7
    figx = ncols*5.4#6.5#4.8
    figy = nrows*3.2#3.1#3.7
    #figx = ncols*6.05#4.8
    #figy = nrows*3.35#3.7
    if single == True:
        ncols = 1
        figx = 6
        figy = 4.5
    if 'samey' in kwargs:
        ymins = []
        ymaxs = []
        for i in range(0,dimx):
            item = dataframe.iloc[-i-1]
            #item = dataframe.iloc[i]
            ymodel = item.model/max(item.model)
            yobs = item.obs/max(item.model)
            ymin = min(min(ymodel),min(yobs))
            ymax = max(max(ymodel),max(yobs))
            ymins.append(ymin)
            ymaxs.append(ymax)
        ymin_glob = min(ymins)
        ymax_glob = max(ymaxs)
        yrange_glob = ymax_glob-ymin_glob
    if 'samex' in kwargs:
        xmin_glob = min([min(dataframe.iloc[i].wave - dataframe.iloc[i].wave_cen) for i in range(0,dimx)])
        xmax_glob = max([max(dataframe.iloc[i].wave - dataframe.iloc[i].wave_cen) for i in range(0,dimx)])
        xrange_glob = xmax_glob-xmin_glob

    fig = plt.figure(figsize=(figx, figy))
    for i in range(0,dimx):
        item = dataframe.iloc[-i-1]
        #item = dataframe.iloc[i]
        label = item.element
        label_new = label
        print(label_new)
        #if label in replace_dict:
            #label_new = label.replace(label, replace_dict[label])
        #else: label_new = label
        ax = fig.add_subplot(nrows,ncols,i+1)
        if 'samey' in kwargs:
            ax.set_ylim(ymin_glob-0.01*yrange_glob, ymax_glob+0.01*yrange_glob)
        #ax.set_ylim(0.50,1.15)

        if 'samex' in kwargs:
            ax.set_xlim(xmin_glob-0.01*xrange_glob, xmax_glob+0.01*xrange_glob)
        #"""
        if len(wave)>0 and len(flux)>0:
            xmin = min(item.wave)
            xmax = max(item.wave)
            xran = xmax - xmin
            print(xmin,xmax)
            mask = np.where((xmin<wave)&(wave<xmax))
            fnorm = norm_lin(wave,flux,xmin,xmin+0.01*xran,xmax-0.01*xran,xmax)
            #kwargs = {"applyto":flux/fnorm}
            fnorm_tospas = norm_lin(item.wave,item.model/max(item.model),
                                    xmin,xmin+0.01*xran,
                                    xmax-0.01*xran,xmax)
            print(len(item.wave))
            print(len(wave))
            print(len(fnorm_tospas))
            fnorm_tospas = np.interp(wave, item.wave, fnorm_tospas)
            # fnorm = 1/max(flux_cut)
            wave_cut = wave[mask]
            flux_norm = flux/fnorm*fnorm_tospas
            #flux_norm = flux/fnorm
            flux_plot  = flux_norm[mask]
            ax.plot(wave_cut - item.wave_cen,
                    flux_plot,
                    lw=0.5,color='blue')
        #"""
        ax.plot(item.wave - item.wave_cen,
                 item.obs/max(item.model),
                 lw=0.5,color='black', zorder=1)
        ax.plot(item.wave - item.wave_cen,
                 item.model/max(item.model),
                 lw=0.85,color='red',
                 label=label_new)
        plt.annotate(s=label_new,size=13,xy=(0, 1),
                     xycoords='axes fraction',
                     xytext=(0.02, 0.015),
                     textcoords='axes fraction',
                     ha='left',va='bottom',color='black')
        plt.annotate(s=r'$\lambda _0$'+str(int(item.wave_cen)),
                     size=13,xy=(0, 1),
                     xycoords='axes fraction',
                     xytext=(0.98, 0.015),
                     textcoords='axes fraction',
                     ha='right',va='bottom',color='black')
        ax.tick_params(axis='both',which='major',direction='in',
                       top=True, bottom=True, left=True, right=True)

    if single == True:
        plt.xlabel("$\lambda -\lambda _0$ ($\AA$)")
        plt.ylabel(r"Normalized Flux")
        plt.tight_layout()
    else:
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none',
                    top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.xlabel("$\lambda -\lambda _0$ ($\AA$)",fontsize=13, labelpad=15)
        plt.ylabel(r"Normalized Flux",fontsize=13, labelpad=13)
        plt.tight_layout(pad=0,w_pad=1,h_pad=0.2)

    plt.close()
    return fig


def main():
    import os, argparse
    import numpy as np
    import pandas as pd
    from matplotlib.backends.backend_pdf import PdfPages
    # from plot_normpanels_resid.py import norm_lin
    
    description = '''\
    description:
      plot SPAS output files
        '''
    epilog = '''\
    examples:
      plot_outspas.py DAO_FUSE_WD1214+267_3D.out
      plot_outspas.py DAO_FUSE_WD1214+267_3D.out --single --save_space --ruler
      plot_outspas.py DAO_FUSE_WD1214+267_3D.out -cs synspec.7 48000
      plot_outspas.py 3D_atm.out -ss -r -ng -sy -rl 0.2
        '''
    parser = argparse.ArgumentParser(
              prog='plot_outspas.py',
              formatter_class=argparse.RawDescriptionHelpFormatter,
              description=description,
              epilog=epilog)
    parser.add_argument("-s","--single",
                        action='store_true', default=False,
                        help="Also save single subplots")
    parser.add_argument("-ss","--save_space",
                        action='store_true', default=False,
                        help="Minimize space between lines")
    parser.add_argument("-r","--ruler",
                        action='store_true', default=False,
                        help="Plot ruler")
    parser.add_argument("-rl","--ruler_length",
                        nargs='?', default=False,type=float,
                        help="Fix ruler length")
    parser.add_argument("-cs","--comp_syn",
                        nargs=2, default=False,
                        help="Compare with add. synthetic model; model R vrad")
    parser.add_argument("-ng","--nogrid",
                        action='store_true', default=False,
                        help="Don't plot grid lines")
    parser.add_argument("-sa","--sortatm",
                        action='store_true', default=False,
                        help="Seperate atm. plots by H/HeII and HeI")
    parser.add_argument("-sy","--samey",
                        action='store_true', default=False,
                        help="Same y-axis in mosaik plot")
    parser.add_argument("-sx","--samex",
                        action='store_true', default=False,
                        help="Same y-axis in mosaik plot")
    parser.add_argument("-ro","--reorder",
                        action='store_true', default=False,
                        help="Indstead of sorting windows by central wavelength,\
                        sort by window names (e.g. Ha before Hb or 1_Hb before 2_Hb)")
    parser.add_argument("-me","--mask_element",
                        nargs='?', default=False,
                        help="Only show windows that contain string in name")
    parser.add_argument("-mr","--mask_range",
                        nargs=2, default=False,type=float,
                        help="Only show windows with mask_range[0] < wave_cen < mask_range[1]")
    parser.add_argument("-c","--center",
                        action='store_true', default=False,
                        help="Align plots at central wavelength (not SPAS wavelength),\
                        sort by window names (e.g. Ha before Hb or 1_Hb before 2_Hb)")
    parser.add_argument('-sb',"--smooth_box",
                        nargs=1, default=[0],type=float,
                        help="Smooth spectra using box filter of specified width in pixels")
    parser.add_argument('argv', nargs='*')
    args = parser.parse_args()
    print(args)

    spasout = str(args.argv[0])
    filename =  os.path.basename(spasout)

    nrows = rawincount(spasout)
    par,df = read_rpar(nrows,spasout)

    with open(spasout, 'r') as file:
        data    = file.readlines()

    waves,flux_obss,flux_models = read_ranges(data,par)

    if args.smooth_box[0] >= 2:
        from scripts_instrument import (est_SNR,
                               convolve_box)
        npix = int(args.smooth_box[0])
        print("applying boxcar smoothing with",npix,"pixels")
        flux_obss = [convolve_box(item,npix)[2:-2] for item in flux_obss]
        flux_models = [convolve_box(item,npix)[2:-2] for item in flux_models]
        waves = [item[2:-2] for item in waves]
        #error_obs, SNR = est_SNR_array(wave_obs, flux_obs, 1)

    df = df.assign(wave=pd.Series(waves),
                       obs=pd.Series(flux_obss),
                       model=pd.Series(flux_models))

    for index, row in df.iterrows():
        print(row['element'],min(row['wave']), max(row['wave']))

    print('len(waves)',len(waves))
    if args.mask_range:
        mask_range = (args.mask_range[0] < df.wave_cen) &\
                     ( df.wave_cen < args.mask_range[1])
        df = df[mask_range]

    if args.mask_element:
        df = df[df.element.str.contains(args.mask_element)]
    if args.reorder:
        # sort values by string
        # remove numbers at beginning of string
        newstrings = []
        numbers = []
        print("df['element']\n",df['element'])
        for item in df['element']:
            print(item)
            number = item.split('_',1)[0]
            newstr = item.split('_',1)[-1]
            try:
                number = int(number)
                print(number)
            except ValueError:
                print("Warning: Ordering by string (not int)")
            print(number,newstr)
            newstrings.append(newstr)
            numbers.append(number)
        df['element'] = newstrings
        df['numbers'] = numbers
        print("df['element']\n",df['element'])
        print("df['numbers']\n",df['numbers'])
        df = df.sort_values(by=['numbers'], ascending=False).reset_index(drop=True)
    else:
        df = df.sort_values(by=['wave_cen'], ascending=False).reset_index(drop=True)
        if args.sortatm:
            df = df.sort_values(by=['wave_cen'], ascending=True).reset_index(drop=True)

    data_edit = [x.replace('# ','') for x in data[2:8]]
    data_string_1 = ''.join(data_edit[:3])
    data_string_2 = ''.join(data_edit[3:])

    if args.comp_syn:
        synname = args.comp_syn[0]
        resolution = float(args.comp_syn[1])
        wave_syn, flux_syn = np.loadtxt(synname,
                                        usecols=(0,1),
                                        unpack=True,
                                        dtype=float)
        c = 299792.458#km/s
        v_rad = float(args.comp_syn[2])
        fac_toobs = np.sqrt((1+(v_rad/c))/(1-(v_rad/c)))
        wave_syn = np.multiply(wave_syn,fac_toobs)
        if synname.startswith("R") == False:
            wave_syn,flux_syn = edit_sec(wave_syn,flux_syn,resolution)
        print(wave_syn[:10])
        print(flux_syn[:10])

    plot_multicolumn_keys, plot_df_keys, plot_mosaic_keys = {}, {}, {}
    if args.sortatm:
        plot_multicolumn_keys['sortatm'] = True
    if args.center:
        plot_multicolumn_keys['center'] = True
        plot_df_keys['center'] = True
    if args.nogrid:
        plot_multicolumn_keys['nogrid'] = True
        plot_df_keys['nogrid'] = True
    if args.ruler:
        plot_multicolumn_keys['ruler'] = True
        plot_df_keys['ruler'] = True
        if args.ruler_length:
            plot_multicolumn_keys['ruler_length'] = args.ruler_length
            plot_df_keys['ruler_length'] = args.ruler_length
    if args.save_space:
        plot_multicolumn_keys['save_space'] = True
        plot_df_keys['save_space'] = True
    if args.comp_syn:
        plot_multicolumn_keys.update({"wave":wave_syn, "flux":flux_syn})
    if args.samey:
        plot_mosaic_keys['samey'] = True
    if args.samex:
        plot_mosaic_keys['samex'] = True

    print("plot_singlemosaic()")
    plot1 = plot_singlemosaic(df,2,**plot_mosaic_keys)
    #print("plot_singlecolumn()")
    #plot0 = plot_singlecolumn(df,**plot_df_keys)
    #print("plot_multicolumn()")
    #plot2 = plot_multicolumn(df,10,**plot_multicolumn_keys)

    try:
        os.makedirs('pdf/')
    except FileExistsError:
        print('directory exists:','pdf/')

    savename=['pdf/'+filename+'.pdf','']

    pp = PdfPages(savename[0])
    plotm = meta_data(42070,5.612,0.332,-17,2,0,data_edit)
    pp.savefig(plotm)
    pp.savefig(plot1)
    #pp.savefig(plot0)
    #pp.savefig(plot2)
    if args.single:
        kwargs = {"single":True}
        if args.comp_syn:
            kwargs.update({"wave":wave_syn, "flux":flux_syn})
        for i in range(0,len(df)):
            plot_single = plot_singlemosaic(df.iloc[[i]],1,**kwargs)
            pp.savefig(plot_single)
    pp.close()
    print('# fit windows:',len(waves))
    print('saved as:',savename[0])

if __name__ == "__main__":
    main()
