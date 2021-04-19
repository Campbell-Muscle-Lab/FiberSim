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

    results_file_strings = ['../sim_output/results_0.txt',
                            '../sim_output/results_1.txt',
                            '../sim_output/results_2.txt']
    output_file_string = 'c:/temp/stretches.png'

    # Make a figure
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches([9, 9])
    spec = gridspec.GridSpec(nrows=6, ncols=3, figure=fig,
                             wspace=1)
    axs=[]
    for r in range(0,6):
        for c in range(0,3):
            axs.append(fig.add_subplot(spec[r,c]))

    d0 = pd.read_csv(results_file_strings[0], sep='\t')
    
    for c in range(0,3):
        d = pd.read_csv(results_file_strings[c], sep='\t')
        x = d['time']

        axs[0+(c)].plot(x, d['pCa'])
        axs[0+(c)].set_ylim([9.5, 4])
        axs[0+(c)].set_ylabel('pCa')

        axs[3+(c)].plot(x, d['hs_length'])
        axs[3+(c)].set_ylim([1150, 1350])
        axs[3+c].set_ylabel('HS length')
        
        axs[6+c].plot(x, d['force'], label='Total')
        axs[6+c].plot(x, d['titin_force'], label='Titin')
        axs[6+c].set_ylim([0, 200000])
        axs[6+c].set_ylabel('Force')
        if (c==2):
            axs[6+c].legend(loc='upper left',
                             bbox_to_anchor=[1.05, 1])
        
        axs[9+c].plot(x, d['a_pop_0'], label='Off')
        axs[9+c].plot(x, d['a_pop_1'], label='On')
        axs[9+c].set_ylim([0, 1])
        axs[9+c].set_ylabel('Thin')
        if (c==2):
            axs[9+c].legend(loc='upper left',
                             bbox_to_anchor=[1.05, 1])
        
        
        axs[12+c].plot(x, d['m_pop_0'], label='SRX')
        axs[12+c].plot(x, d['m_pop_1'], label='DRX')
        axs[12+c].plot(x, d['m_pop_2'], label='Pre-power')
        axs[12+c].plot(x, d['m_pop_3'], label='Post-power')
        axs[12+c].set_ylim([0, 1])
        axs[9+c].set_ylabel('Thick')
        if (c==2):
            axs[12+c].legend(loc='upper left',
                             bbox_to_anchor=[1.05, 1])
        
        axs[15+c].plot(x, d['c_pop_0'], label='Detached')
        axs[15+c].plot(x, d['c_pop_1'], label='Attached')
        axs[15+c].set_ylim([0, 1])
        axs[15+c].set_ylabel('MyBPC')
        if (c==2):
            axs[15+c].legend(loc='upper left',
                             bbox_to_anchor=[1.05, 1])

    for c in range(0,3):
        d = d0
        x = d['time']

        axs[0+(c)].plot(x, d['pCa'])
        axs[0+(c)].set_ylim([9.5, 4])
        axs[0+(c)].set_ylabel('pCa')

        axs[3+(c)].plot(x, d['hs_length'], 'k-')
        
        axs[6+c].plot(x, d['force'], 'k-')


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
