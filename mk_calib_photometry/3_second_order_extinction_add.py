#! /usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
####          Configuration: modify the file in this section            ####
############################################################################

filter_list = [
    'B',
    'V',
    ]

table_list = [
    'output_GD_279_add/tables/trans_para_GD_279.dat',
    'output_NGC188_add/tables/trans_para_NGC_188.dat',
    'output_NGC752_add/tables/trans_para_NGC_752.dat',
    'output_NGC884_add/tables/trans_para_NGC_884.dat',
    'output_NGC6939_add/tables/trans_para_NGC_6939.dat',
    'output_NGC6940_add/tables/trans_para_NGC6940.dat',
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

outdir = 'output_add_all'

weights = False

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
