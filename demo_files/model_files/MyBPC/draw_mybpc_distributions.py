# -*- coding: utf-8 -*-
"""
Created on Sat May  7 16:02:27 2022

@author: kscamp3
"""

import os
import json

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def draw_mybpc_distributions():
    """ function draws the number of mybpc pointing to each thin filament
        for models with different c_sd_angle_deviations """
    
    # Variables
    
    base_status_folder = 'sim_output/hs';
    condition_folders = ['model_c_sd_angle_deviation_absent',
                         'model_c_sd_angle_deviation_0',
                         'model_c_sd_angle_deviation_30']
    status_file = 'hs_1_time_step_1.json'
    output_file_string = 'sim_output/mybpc_histogram.png'
    
    # Make a folderrue
    fig = plt.figure(constrained_layout=False)
    fig.set_size_inches([6,6])
    fig_rows = 3
    spec=gridspec.GridSpec(nrows=fig_rows, ncols=1, figure=fig,
                           left=0.3, hspace=1)
    
    ax=[]
    for i in range(fig_rows):
        ax.append(fig.add_subplot(spec[i,0]))
    
    # Loop through the condition folders
    for i, cond_fold in enumerate(condition_folders):
        
        # Generate the filename
        file_name = os.path.join(base_status_folder,
                                 cond_fold,
                                 status_file)
        
        # Load as dict
        with open(file_name, 'r') as f:
            st = json.load(f)
            
        # Loop through the thick filaments, collating the nearest a_f values
        nearest_a_f = []
        thick = st['thick']
        
        for th in thick:
            nearest_a_f = nearest_a_f + th['pc_nearest_a_f']
            
        # Determine the number of thin filaments
        no_of_thin_filaments = len(st['thin'])
        
        x = list(range(1,no_of_thin_filaments+1))
       
        # Calculate the histrogram
        [counts, edges]= np.histogram(nearest_a_f, x)
        
        # Plot
        ax[i].bar(edges[0:-1], counts, facecolor='b')
        ax[i].set_ylabel('Number of\nMYBPC nearest\nthin filament',
                         rotation=0,
                         verticalalignment='center',
                         labelpad=40)
        ax[i].set_title(cond_fold)
        
        if i==2:
            ax[i].set_xlabel('Thin filament number')
            
    print('Saving figure to: %s' % output_file_string)
    fig.savefig(output_file_string)
            
    
if __name__ == "__main__":
    draw_mybpc_distributions()
