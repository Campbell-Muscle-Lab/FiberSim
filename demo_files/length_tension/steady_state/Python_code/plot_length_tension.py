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

def plot_length_tension():
    
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
  
    # Now we can make a figure
    no_of_rows = 1
    no_of_cols = 1
    
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([5, 7])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
    
    # Arrange data holders
    hsl = []
    pas_force = []
    act_force = []
    tot_force = []
            
    # Loop through the simulations
    for i in range(no_of_conditions):
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % (i+1)))
        
        print('\nCondition_folder: %s\n' % condition_folder)

        # Get the sim_data files for the condition
        for file in os.listdir(condition_folder):
            
            if (file.startswith(('rates','status'))):
                continue
            else:
                # Adjust path
                file = os.path.join(condition_folder, file)

            # Open up the file                    
            d = pd.read_csv(file, sep='\t')
            
            # Calculate some numbers
            hs_length = d['hs_1_length'].iloc[0]
            passive_force = d['hs_1_force'].iloc[0]
            total_force = np.mean(d['hs_1_force'].iloc[-10:-1])
            active_force = total_force - passive_force
            
            # Store
            hsl.append(hs_length)
            pas_force.append(passive_force)
            act_force.append(active_force)
            tot_force.append(total_force)

    # Now plot
    ax[0].plot(hsl, pas_force, 'bo', label='Passive')
    ax[0].plot(hsl, act_force, 'rs', label='Active')
    ax[0].plot(hsl, tot_force, 'gd', label='Total')
    ax[0].plot([hsl[0], hsl[-1]], [0, 0], 'k-')
    ax[0].legend()
    
    # Save fig
    output_file_string = os.path.join(top_data_folder, 'length_tension.png')
    fig.savefig(output_file_string, bbox_inches='tight')
    
if __name__ == "__main__":
    plot_length_tension()
