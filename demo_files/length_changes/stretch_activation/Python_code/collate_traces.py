# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 17:21:31 2025

@author: Campbell
"""

import os
import json
import copy
import re

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def collate_sim_files(pCa_value, folder_indices):
    
    top_data_folder = '../sim_data/sim_output'
    
    sim_files = []
    pCa_string = 'pCa_%i' % (10*pCa_value)
    
    for fi in folder_indices:
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % fi))
        
        all_files = os.listdir(condition_folder)
        
        for af in all_files:
            if (pCa_string in af):
                sim_files.append(os.path.join(condition_folder,
                                              af))
    return (sim_files)
                    

def create_axes():
    # Make a figure
    no_of_rows = 5
    no_of_cols = 1
    
    pCa_row = 1
    hs_force_row = 2
    hs_a_row = 3
    hs_m_row = 4
    
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([no_of_cols * 3, 10])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
            
    return fig, ax

def fig_collate_traces(sim_files, fig, ax):
    
    for sfs in sim_files:
        
        d = pd.read_csv(sfs, sep='\t')
        
        ax[0].plot(d['time'], d['hs_1_pCa'], '-')
        ax[1].plot(d['time'], d['hs_1_force'], '-')
        ax[2].plot(d['time'], d['hs_1_length'], '-')
        
    ofs = 'd:/temp/k2.png'
    print('Saving summary_figure to: %s' % ofs)
    fig.savefig(ofs, bbox_inches='tight')
        
        
        
        
if __name__ == '__main__':

    fig, ax = create_axes();    
    
    sim_files = collate_sim_files(6.3, [1, 2, 3, 4, 5])
    
    print(sim_files)
    
    fig_collate_traces(sim_files, fig, ax)
        
        
        
    
    
    
    