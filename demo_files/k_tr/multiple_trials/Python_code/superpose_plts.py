# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 15:14:27 2021

@author: kscamp3
"""

import os
import json

import numpy as np
import pandas as pd

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

def analyze_simulations():

    results_file_string_base = '../sim_output/results_'
    output_file_string = '../temp/superposed.png'

    # Make a figure
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches([6, 6])
    spec = gridspec.GridSpec(nrows=3, ncols=1, figure=fig,
                             wspace=1)
    axs=[]
    for r in range(0,3):
        for c in range(0,1):
            axs.append(fig.add_subplot(spec[r,c]))

    for i in range(0,4):
        rfs = ('%s%i.txt' % (results_file_string_base, i))

        d = pd.read_csv(rfs, sep='\t')

        x = d['time']

        axs[0].plot(x, d['pCa'])
        axs[0].set_ylim([9.5, 4])
        axs[0].set_ylabel('pCa')

        axs[1].plot(x, d['hs_length'])
        axs[1].set_ylim([1140, 1210])
        axs[1].set_ylabel('HS length')
        
        axs[2].plot(x, d['force'])
        axs[2].set_ylim([0, 120000])


    # Save figure
    print('Saving current fit to %s' % output_file_string)
    # Check folder exists and make it if not
    dir_name = os.path.dirname(os.path.abspath(
        output_file_string))
    if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
    fig.savefig(output_file_string)
    plt.close()


    # # Add legends
    # for i in range(2,4):
    #     axs[i].legend(legend_symbols_N, legend_strings_N,
    #               loc='upper right',
    #               bbox_to_anchor=(1.5, 1))
    # for i in range(4,6):
    #     axs[i].legend(legend_symbols_m, legend_strings_m,
    #               loc='upper right',
    #               bbox_to_anchor=(1.5, 1))


if __name__ == "__main__":
    analyze_simulations()
