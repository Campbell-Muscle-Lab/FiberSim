# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 21:56:54 2024

@author: ken
"""

import os
import sys

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def return_fit(thread_id_dir):
    """ returns a single value defining the least squares fit
        between the current simulation and the target data """
       
    # Variables
    sim_data_file = 'sim_data/sim_output/1/sim_prot_1_1_r1.txt'
    target_data_file = '../../target/target_sim.txt'
    trial_errors_file = 'working/trial_errors.xlsx'
    output_image_folder = 'sim_data/sim_output'
    
    # Code    

    # Adapt because files are relative to the thread_id_path
    
    # Set the trial errors file
    trial_errors_file = os.path.join(thread_id_dir, trial_errors_file)
    
    sim_data_file = Path(os.path.join(thread_id_dir, sim_data_file)).resolve()

    # Open it
    sim_data = pd.read_csv(sim_data_file, sep='\t')
    
    # Open the target_data
    target_data_file = os.path.join(thread_id_dir, target_data_file)
    
    target_data = pd.read_csv(target_data_file, sep='\t')

    # Calculate the error
    dy = (sim_data['hs_1_force'].to_numpy() -
         target_data['hs_1_force'].to_numpy())
    dy_sq = np.power(dy, 2)
    e = np.sum(dy_sq)

    # Write error data to file
    d = dict()
    d['error_cpt_1'] = e
    d['error_total'] = e
    
    # Make a dataframe    
    df = pd.DataFrame(data=d, index=[0])
    print(df)
    
    # Check the dir exists
    worker_parent_dir = Path(trial_errors_file).parent
    if not os.path.isdir(worker_parent_dir):
        os.makedirs(worker_parent_dir)
    
    # Clean the file and then write
    if (os.path.exists(trial_errors_file)):
        os.remove(trial_errors_file)
    df.to_excel(trial_errors_file, index=False)
    
    # Make a figure
    fig = plt.figure(constrained_layout=False)
    spec = fig.add_gridspec(nrows=1, ncols=1,
                      left=0.3, right=0.95, wspace=0.1,
                      bottom = 0.2)
    fig.set_size_inches([3.5, 2.5])
    ax = fig.add_subplot(spec[0, 0])
    
    ax.plot(target_data['time'], target_data['hs_1_force'], 'k-')
    ax.plot(sim_data['time'], sim_data['hs_1_force'], 'b-')
    
    ofs = os.path.join(thread_id_dir, output_image_folder,
                       'trace_comparison.png')
    fig.savefig(ofs, dpi=200, bbox_inches='tight')
    
    
if __name__ == "__main__":
    # Thread path gets included with function call as second argument    
    # Code needs this to work out where files are
    return_fit(sys.argv[1])# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 17:08:18 2025

@author: Campbell
"""

