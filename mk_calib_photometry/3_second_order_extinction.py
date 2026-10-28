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

import time
import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
plt.switch_backend('Agg')

from ost_photometry import checks, terminal_output
from ost_photometry.analyze.calibration.second_order_extinction import (
    run_second_order_campaign,
)

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    start_time = time.time()
    terminal_output.print_header("Calibration 3/3: second order extinction")
    checks.check_output_directories(outdir)

    results = run_second_order_campaign(
        table_list,
        filter_list,
        outdir,
        apply_weights=weights,
    )

    for fit in results:
        terminal_output.print_section(f"Filter {fit.filter}")
        terminal_output.print_to_terminal(
            f"Transformation coefficient  T{fit.legacy_t_label} = "
            f"{fit.intercept_t:.5} +/- {fit.intercept_t_err:.5}",
            indent=2, style_name="OKGREEN")
        terminal_output.print_to_terminal(
            f"Second order extinction     k\"{fit.legacy_t_label} = "
            f"{fit.k_second_order:.5} +/- {fit.k_second_order_err:.5}",
            indent=2, style_name="OKGREEN")

    terminal_output.print_done(start_time, {"Output": outdir})
