#! /usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
####          Configuration: modify the file in this section            ####
############################################################################

filter_list = [
    #'U',
    'B',
    'V',
    #'R',
    ]
C_name      = {
    #'U':['Cuuv', 'Cvuv'],
    'B':['Cbbv', 'Cvbv'],
    'V':['Cbbv', 'Cvbv'],
    #'R':['Crvr', 'Cvvr'],
    }
C_err_name  = {
    #'U':['Cuuv_err', 'Cvuv_err'],
    'B':['Cbbv_err', 'Cvbv_err'],
    'V':['Cbbv_err', 'Cvbv_err'],
    #'R':['Crvr_err', 'Cvvr_err'],
    }

table_list  = [
    'output_GD_279_add/tables/trans_para_GD_279.dat',
    #'output_Melotte111_add/tables/trans_para_Melotte_111.dat',
    'output_NGC188_add/tables/trans_para_NGC_188.dat',
    'output_NGC752_add/tables/trans_para_NGC_752.dat',
    'output_NGC884_add/tables/trans_para_NGC_884.dat',
    'output_NGC6939_add/tables/trans_para_NGC_6939.dat',
    'output_NGC6940_add/tables/trans_para_NGC_6940.dat',
    'output_NGC7142_add/tables/trans_para_NGC_7142.dat',
    'output_NGC7790_add/tables/trans_para_NGC_7790.dat',
    'output_SA29_add/tables/trans_para_SA_29.dat',
    'output_SA95_add/tables/trans_para_SA_95.dat',
    'output_SA98_add/tables/trans_para_SA_98.dat',
    'output_NGC6939_20220117_add/tables/trans_para_NGC_6939.dat',
    'output_NGC188_20220117_add/tables/trans_para_NGC_188.dat',
    'output_NGC7789_20220117_add/tables/trans_para_NGC_7789.dat',
    'output_NGC457_20220117_add/tables/trans_para_NGC_457.dat',
    'output_M37_add/tables/trans_para_M_37.dat',
    ]

outdir      = 'output_add_all'

#   Apply weights in calculation
weights = False
#weights = True

############################################################################
####                            Libraries                               ####
############################################################################

import sys

import time

import numpy as np

from astropy.table import Table

import matplotlib.pyplot as plt
plt.switch_backend('Agg')

from ost_photometry.analyze.aux import (lin_func,
                             fit_curve,
                             )

from ost_photometry import checks

from ost_photometry.style import bcolors

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    ###
    #   Check output directories
    #
    checks.check_out(
        outdir,
        )

    for j, filt in enumerate(filter_list):
        for key, C_list in C_name.items():
            #   Find the entry in the C_name dictionary that belongs to the
            #   current filter
            if key == filt:
                C_err_list = C_err_name[key]
                #   Loop over the C values
                for u, C_entry in enumerate(C_list):
                    C_err_entry = C_err_list[u]

                    #   Define lists for the C values and airmass from the
                    #   different cluster
                    C       = []
                    C_err   = []
                    airmass = []
                    name    = []
                    for i, path in enumerate(table_list):
                        tbl = Table.read(path, format='ascii')
                        if C_entry in tbl.colnames:
                            if 'airmass_'+filt in tbl.colnames:
                                C.append(tbl[C_entry].data[0])
                                C_err.append(tbl[C_err_entry].data[0])
                                airmass.append(tbl['airmass_'+filt].data[0])
                                name.append(tbl['name'].data[0])


                    #   Initial guess for the parameters
                    x0    = np.array([0.0, 0.0])

                    #   Set sigma, unsing errors
                    if weights:
                        sigma = np.array(C_err)
                        #sigma = 1./np.sqrt(C_err)
                        #sigma = np.sqrt(C_err)
                    else:
                        sigma = 0.

                    #   Fit
                    fit_func=lin_func
                    T, T_err, k, k_err = fit_curve(
                        fit_func,
                        airmass,
                        C,
                        x0,
                        sigma,
                        )

                    #   Fit data
                    x_lin = np.sort(airmass)
                    y_lin = fit_func(x_lin, T, k)

                    titel  = "C values vs. air mass"
                    plabel = 'slope = '+str(k)+', k = '+str(-k)\
                             +' +/- '+str(k_err)+'; T'+C_entry[1:4]\
                             +' = '+str(T)+' +/- '+str(T_err)
                    xlabel = 'x [air mass]'
                    ylabel = C_entry
                    path   = outdir+'/C_vs_x_'+filt+'_T'+C_entry[1:4]+'.pdf'

                    ##  Make plot
                    fig = plt.figure(figsize=(20,9))
                    #   Set title
                    fig.suptitle(titel, fontsize=30)
                    #   Plot data
                    plt.errorbar(
                        airmass,
                        C,
                        yerr=C_err,
                        #yerr=sigma,
                        color='blue',
                        marker='.',
                        mew=0.0,
                        linestyle='none',
                        )
                    #   Plot labels
                    for i in range(0, len(name)):
                        plt.annotate(
                            name[i],
                            (airmass[i], C[i]),
                            #textcoords='data',
                            xytext=(airmass[i], C[i]),
                            textcoords='offset points',
                            )
                    #   Plot fit
                    plt.plot(
                        x_lin,
                        y_lin,
                        linestyle='-',
                        color='red',
                        linewidth=0.8,
                        label=plabel,
                        )
                    #   Set legend
                    plt.legend(
                        bbox_to_anchor=(0.,1.02,1.0,0.102),
                        loc=3,
                        ncol=4,
                        mode='expand',
                        borderaxespad=0.,
                        )

                    plt.grid(visible=True, axis='y')

                    #   Set x and y axis lable
                    plt.xlabel(xlabel, fontsize=20)
                    plt.ylabel(ylabel, fontsize=20)
                    #   Save plot
                    plt.savefig(path,bbox_inches='tight',format='pdf')
                    plt.close()

                    indent = '      '
                    #   Print results
                    print(indent+'########################################')
                    print(
                        bcolors.BOLD
                        +indent+'Transformation coefficients ('+filt+'):'
                        +bcolors.ENDC
                        )
                    print(
                        indent+'   T'+C_entry[1:4]+' = ',
                        f"{T:.5}",
                        '+/-',
                        f"{T_err:.5}"
                        )
                    print(
                        bcolors.BOLD
                        +indent+'Second order extinction coefficients ('
                        +filt+'):'
                        +bcolors.ENDC
                        )
                    print(
                        indent+'   k"'+C_entry[1:4]+' = ',
                        f"{-k:.5}",
                        '+/-',
                        f"{k_err:.5}"
                        )
                    print(indent+'########################################')


