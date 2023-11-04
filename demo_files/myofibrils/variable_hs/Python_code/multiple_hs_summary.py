# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:16:23 2023

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

from natsort import natsorted

from pathlib import Path

def summary_figure():
    
    top_data_folder = '../sim_data/sim_output'
    
    # Adapt because it is relative to this file
    parent_dir = Path(__file__).parent.absolute()
    top_data_folder = Path(os.path.join(parent_dir, top_data_folder)).resolve()
    
    print(top_data_folder)
    
    # Find out how many simulations there are
    no_of_conditions = 1
    keep_going = True

    while (keep_going):
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % (no_of_conditions)))

        if os.path.isdir(condition_folder):
            no_of_conditions = no_of_conditions + 1
        else:
            keep_going = False
            
    no_of_conditions = no_of_conditions - 1

    # Now we can make a figure with one column per condition
    no_of_rows = 4
    no_of_cols = no_of_conditions
    
    pCa_row = 1
    hs_force_row = 2
    hs_titin_row = 3
    hs_length_row = 4
    
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([no_of_conditions * 3, 10])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
    
    min_force = np.inf
    max_force = -np.inf
    
    min_titin = np.inf
    max_titin = -np.inf
    
    min_hsl = np.inf
    max_hsl = -np.inf
    
    # Get a color map
    cmap = matplotlib.colormaps['viridis']
    
            
    # Loop through the simulations
    for i in range(no_of_conditions):
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % (i+1)))
        
        for file in os.listdir(condition_folder):
            if (file.startswith('rates')):
                continue
            
            file = os.path.join(condition_folder, file)
        
            d = pd.read_csv(file, sep='\t',
                            na_values=['-nan(ind'])
            
            # Now get the number of half-sarcomeres
            no_of_hs = 0
            for h in d.columns.to_list():
                if ('pCa' in h):
                    n = [float(num) for num in re.findall(r'[\d]+', h)]
                    if (n[0] > no_of_hs):
                        no_of_hs = int(n[0])
                        
            # Set indices for colormap
            cmap_values = np.linspace(0, 1, no_of_hs)
                        
            # pCa
            plot_index = ((pCa_row - 1) * no_of_cols) + i
            ax[plot_index].plot(d['time'], d['hs_1_pCa'], '-')
                   
            # hs_force
            plot_index = ((hs_force_row-1) * no_of_cols) + i
            for j in range(no_of_hs):
                col_string = 'hs_%i_force' % (j+1);
                ax[plot_index].plot(d['time'], d[col_string], '-',
                                    color=cmap(cmap_values[j]))
                if (d[col_string].max() > max_force):
                    max_force = d[col_string].max()
                if (d[col_string].min() < min_force):
                    min_force = d[col_string].min()
                    
            # hs_titin
            plot_index = ((hs_titin_row -1) * no_of_cols) + i
            for j in range(no_of_hs):
                col_string = 'hs_%i_titin' % (j+1)
                ax[plot_index].plot(d['time'], d[col_string], '-',
                                    color=cmap(cmap_values[j]))
                if (d[col_string].max() > max_titin):
                    max_titin = d[col_string].max()
                if (d[col_string].min() < min_titin):
                    min_titin = d[col_string].min()
                        
            # hs_length
            plot_index = ((hs_length_row-1) * no_of_cols) + i
            for j in range(no_of_hs):
                col_string = 'hs_%i_length' % (j+1);
                ax[plot_index].plot(d['time'], d[col_string], '-',
                                    color=cmap(cmap_values[j]))
                if (d[col_string].max() > max_hsl):
                    max_hsl = d[col_string].max()
                if (d[col_string].max() < min_hsl):
                    min_hsl = d[col_string].min()
                
    # Apply limits
    for i in range(no_of_conditions):
        
        # pCa
        plot_index = ((pCa_row - 1) * no_of_cols) + i
        ax[plot_index].set_ylim([9.5, 4])
        
        # Force
        plot_index = ((hs_force_row - 1) * no_of_cols) + i
        ax[plot_index].set_ylim([0, max_force + 5000])
        
        # Force
        plot_index = ((hs_titin_row - 1) * no_of_cols) + i
        ax[plot_index].set_ylim([min_titin - 5000, max_titin + 5000])
        
        # HSL
        plot_index = ((hs_length_row - 1) * no_of_cols) + i
        ax[plot_index].set_ylim([min_hsl - 50, max_hsl + 50])
        
            
    ofs = os.path.join(top_data_folder, 'summary.png')
    print('Saving summary_figure to: %s' % ofs)
    fig.savefig(ofs, bbox_inches='tight')
            

if __name__ == "__main__":
    summary_figure()
    