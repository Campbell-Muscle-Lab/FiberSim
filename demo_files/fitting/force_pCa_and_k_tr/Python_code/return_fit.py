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
            
        if (i==1):
            test_value = pCa_data['y_min'].iloc[0]
            target_value = target_data['Target'].iloc[1]
            w = target_data['Weight'].iloc[1]
            
        if (i==2):
            test_value = pCa_data['pCa_50'].iloc[0]
            
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

        # Calculate error_component    
        error_components[i] = w * np.power(
            (target_value - test_value) / target_value, 2.0)

    # Calculate final value
    e = np.sum(error_components)    
    
    # Display
    print('Error components')
    print(error_components)
            
    print('Error')
    print(e)

    
if __name__ == "__main__":
    return_fit()