# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 17:48:35 2024

@author: Campbell
"""

"""
Created on Sat Sep 16 15:16:23 2023

@author: Campbell
"""
import os
import sys

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from natsort import natsorted

from pathlib import Path

# Add FiberSim code to path
sys.path.append('../../../code/fiberpy/fiberpy/package/modules/analysis')
import dump_file_analysis as dump_analysis

def plot_myosin_isotypes():
    
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
    no_of_rows = 5
    no_of_cols = no_of_conditions
    
    pCa_row = 1
    hs_force_row = 2
    hs_length_row = 3
    
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
    
    min_hsl = np.inf
    max_hsl = -np.inf
    
            
    # Loop through the simulations
    for i in range(no_of_conditions):
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % (i+1)))
        
        print('\nCondition_folder: %s\n' % condition_folder)

        # Get the sim_data files for the condition
        sim_data_files = []
        status_folders = []
        for file in os.listdir(condition_folder):
            if (file.startswith('rates')):
                continue
            else:
                # Adjust path
                file = os.path.join(condition_folder, file)
                if (os.path.isfile(file)):
                    sim_data_files.append(file)
                else:
                    status_folders.append(file)
                
        # Sort the files in pCa order
        sim_data_files = natsorted(sim_data_files, reverse=True)
        
        # Get the default colors
        color_map = [p['color'] for p in plt.rcParams['axes.prop_cycle']]

        # Now open up the sim_data_file
        for (j,sdf) in enumerate(sim_data_files):
            d = pd.read_csv(sdf, sep='\t')
            
            # Plot pCa
            plot_index = i + ((pCa_row-1) * no_of_cols)
            ax[plot_index].plot(d['time'], d['hs_1_pCa'], '-',
                                color = color_map[i])
            
            # Plot force
            plot_index = i + ((hs_force_row-1) * no_of_cols)
            ax[plot_index].plot(d['time'], d['hs_1_force'], '-',
                                color = color_map[i])
            
            # Now work out the dump files, there is 1 for every time-point
            dump_files = []
            for file in os.listdir(status_folders[i]):
                dump_files.append(file)
            # Sort them
            dump_files = natsorted(dump_files)
            
            # Loop through them
            for (k, file) in enumerate(dump_files):
                file = os.path.join(status_folders[i], file)
                dump_analysis.extract_dump_data(file)
                
                break
    
    # Save fig
    output_file_string = os.path.join(top_data_folder, 'myosin_isotypes.png')
    fig.savefig(output_file_string, bbox_inches='tight')
    
    a = dump_analysis.parse_dump_file_string('ken')
    
    print(a)
    
if __name__ == "__main__":
    plot_myosin_isotypes()
