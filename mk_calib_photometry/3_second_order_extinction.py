#! /usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
####          Configuration: modify the file in this section            ####
############################################################################

filter_list = [
    'U',
    'B',
    'V',
    'R',
    ]

# Paths to per-field tables (``trans_para_*.dat`` or ``trans_para_*.json``)
table_list = [
    'output_GD_279/final/tables/trans_para_GD_279.dat',
    'output_Melotte111_1/final/tables/trans_para_Melotte111.dat',
    'output_NGC7790/final/tables/trans_para_NGC7790.dat',
    'output_SA29/final/tables/trans_para_SA29.dat',
    'output_SA95/final/tables/trans_para_SA95.dat',
    'output_SA98/final/tables/trans_para_SA98.dat',
    'output_NGC188/final/tables/trans_para_NGC_188.dat',
    'output_NGC7142/final/tables/trans_para_NGC_7142.dat',
    'output_NGC6940/final/tables/trans_para_NGC6940.dat',
    'output_NGC6939/final/tables/trans_para_NGC6939.dat',
    'output_NGC752/final/tables/trans_para_NGC752.dat',
    ]

outdir = 'output_all'

weights = True

############################################################################
####                            Libraries                               ####
############################################################################

import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
plt.switch_backend('Agg')

from ost_photometry import checks
from ost_photometry.analyze.calibration.second_order_extinction import (
    run_second_order_campaign,
)
from ost_photometry.style import Bcolors as bcolors

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    checks.check_output_directories(outdir)

    results = run_second_order_campaign(
        table_list,
        filter_list,
        outdir,
        apply_weights=weights,
    )

    indent = '      '
    for fit in results:
        print(indent + '########################################')
        print(
            bcolors.BOLD
            + indent + 'Transformation coefficients (' + fit.filter + '):'
            + bcolors.ENDC
        )
        print(
            indent + '   T' + fit.legacy_t_label + ' = ',
            f"{fit.intercept_t:.5}",
            '+/-',
            f"{fit.intercept_t_err:.5}",
        )
        print(
            bcolors.BOLD
            + indent + 'Second order extinction coefficients ('
            + fit.filter + '):'
            + bcolors.ENDC
        )
        print(
            indent + '   k"' + fit.legacy_t_label + ' = ',
            f"{fit.k_second_order:.5}",
            '+/-',
            f"{fit.k_second_order_err:.5}",
        )
        print(indent + '########################################')
