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

def return_fit():
    """ returns a single value defining the least squares fit
        between the current simulation and the target data """
        
    # Variables
    top_data_folder = '../sim_data/sim_output'
    pCa_data_file = 'pCa_analysis.xlsx'
    pCa_sheet = 'curve_1'
    k_tr_data_file = 'k_tr_analysis.xlsx'
    target_data_file = '../target/target_data.xlsx'
    trial_errors_file = '../working/trial_errors.xlsx'
    comparison_file = '../sim_data/sim_output/comparison.xlsx'
    
    # Code    
    # Adapt because files are relative to this file
    parent_dir = Path(__file__).parent.absolute()
    top_data_folder = Path(os.path.join(parent_dir, top_data_folder)).resolve()
    
    # Open the pCa data
    pCa_data_file = os.path.join(top_data_folder, pCa_data_file)
    
    pCa_data = pd.read_excel(pCa_data_file,
                             sheet_name = pCa_sheet)
    
    # Open the k_tr file
    k_tr_data_file = os.path.join(top_data_folder, k_tr_data_file)
    
    k_tr_data = pd.read_excel(k_tr_data_file)
    
    # Open the target_data
    target_data_file = os.path.join(parent_dir, target_data_file)
    
    target_data = pd.read_excel(target_data_file)
    
    # Set the trial errors file
    trial_errors_file = os.path.join(parent_dir, trial_errors_file)
    
    # Set the comparison file
    comparison_file = os.path.join(parent_dir, comparison_file)
    
    # Display
    print(pCa_data)
    print(k_tr_data)
    print(target_data)
    
    # Build up the error signal
    no_of_error_components = 7
    error_components = np.NaN * np.ones(no_of_error_components)
    
    # Cycle through them
    for i in range(no_of_error_components):
        
        # Pull test values from simulation spreadsheets
        if (i==0):
            test_value = pCa_data['y_min'].iloc[0] + \
                pCa_data['y_amp'].iloc[0]
            
            # Add in a column for test values
            target_data['Test_value'] = np.NaN * np.ones(len(target_data['Target']))
            
        if (i==1):
            test_value = pCa_data['y_min'].iloc[0]
            target_value = target_data['Target'].iloc[1]
            w = target_data['Weight'].iloc[1]
            
        if (i==2):
            test_pCa = pCa_data['pCa_50'].iloc[0]
            test_value = np.power(10, -test_pCa)
            
        if (i==3):
            test_value = pCa_data['n_H'].iloc[0]
            
        if (i==4):
            test_value = k_tr_data['k_tr'].iloc[0]
            
        if (i==5):
            test_value = k_tr_data['k_tr'].iloc[1]

        if (i==6):
            test_value = k_tr_data['k_tr'].iloc[2]

        # Pull targets
        target_value = target_data['Target'].iloc[i]
        w = target_data['Weight'].iloc[i]
        
        # Add test_value to target
        target_data.loc[i, 'Test_value'] = test_value

        # Calculate error_component    
        error_components[i] = w * np.power(
            (target_value - test_value) / target_value, 2.0)

    # Calculate final value
    e = np.sum(error_components)    
    
    # Write error data to file
    d = dict()
    for i in range(len(error_components)):
        d['error_cpt_%i' % (i+1)] = error_components[i]
    d['error_total'] = e
    
    # Make a dataframe    
    df = pd.DataFrame(data=d, index=[0])
    
    print(target_data)
    print(df)
    
    # Check the dir exists
    worker_parent_dir = Path(trial_errors_file).parent
    if not os.path.isdir(worker_parent_dir):
        os.makedirs(worker_parent_dir)
    
    # Clean the file and then write
    if (os.path.exists(trial_errors_file)):
        os.remove(trial_errors_file)
    df.to_excel(trial_errors_file, index=False)
        
    # Now write the test_values file
    target_data.to_excel(comparison_file, index=False)
    
if __name__ == "__main__":
    return_fit()