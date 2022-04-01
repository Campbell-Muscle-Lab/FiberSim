# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 18:38:37 2022

@author: kscamp3
"""

import os

import numpy as np
import pandas as pd

def write_protocol_to_file(prot, prot_file_string):
    """ Writes a protocol defined as a Pandas dataframe to a file_string
        defined with an absolute path """
    
    # Parse the file_string, and check the parent folder exists
    # If not, create the directory
    parent_dir = os.path.dirname(prot_file_string)
    
    if not os.path.isdir(parent_dir):
        print('Creating parent dir: %s' % parent_dir)
        os.makedirs(parent_dir)
    
    # Write prot to file
    print('Writing protocol to: %s' % prot_file_string)
    prot.to_csv(prot_file_string, index=False, sep='\t')

def create_length_control_protocol(time_step=0.0001, n_points=1000,
                              initial_pCa=9.0, step_pCa=4.5,
                              step_pCa_s=0.10,
                              mode_vector = [],
                              delta_hsl=[]):
    """ Creates a length control protocol as a Pandas dataframe """
    
    dt = time_step * np.ones(n_points)
    t = np.cumsum(dt)
    pCa = initial_pCa * np.ones(n_points)
    pCa[t > step_pCa_s] = step_pCa
    
    # Set dhsl as 0 if undefined (thus isometric)
    # Otherwise, fill based on delta_hsl
    if (delta_hsl == []):
        dhsl = np.zeros(n_points)
    else:
        dhsl = delta_hsl
    
    # Set mode_vector as -2 if undefined (thus length control)
    if (mode_vector == []):
        mode = -2 * np.ones(n_points)
    else:
        mode = mode_vector
           
    
    # Assemble
    d = {'dt': dt, 'pCa': pCa, 'dhsl': dhsl, 'mode': mode}
    df = pd.DataFrame(data=d)
    
    return df

def create_force_control_protocol(time_step=0.0001, n_points=1000,
                              initial_pCa=9.0, step_pCa=4.5,
                              step_pCa_s=0.01,
                              mode_vector=[],
                              iso_f=0,
                              iso_start_s=0.08):
    """ Creates a force control protocol as a Pandas dataframe """
    
    dt = time_step * np.ones(n_points)
    t = np.cumsum(dt)
    pCa = initial_pCa * np.ones(n_points)
    pCa[t > step_pCa_s] = step_pCa
    dhsl = np.zeros(n_points)
    if (mode_vector == []):
        mode = -2 * np.ones(n_points)
        mode[t > iso_start_s] = iso_f
    else:
        mode = mode_vector
    
    # Assemble
    d = {'dt': dt, 'pCa': pCa, 'dhsl': dhsl, 'mode': mode}
    df = pd.DataFrame(data=d)
    
    return df

    
    
    
    
    
    
    


