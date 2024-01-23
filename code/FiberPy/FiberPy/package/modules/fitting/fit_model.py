# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:36:17 2024

@author: Campbell
"""

import os
import json
import shutil
import copy

import numpy as np
import pandas as pd

from pathlib import Path

from ..batch import batch

def fit_model(json_analysis_file_string):
    """ Code takes a setup fiel with a model/fitting section, and
        tries to fit the model to data """
    
    # Check the analysis file
    if (not json_analysis_file_string):
        print('fitting: fit_model: no analysis file specified')
        exit(1)
        
    # Load it
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        fitting_struct = json_data['FiberSim_setup']['model']['fitting']

    # Deduce the number of parameters
    no_of_parameters = len(fitting_struct['adjustments'])
    
    if ('initial_guess' in fitting_struct):
        p = np.asarray(fitting_struct['initial_guess'])
    else:
        p = 0.5 * np.ones(no_of_parameters)
        
    # Call the worker
    worker(p, json_analysis_file_string)
    
def worker(p_vector, json_analysis_file_string):
    """ Code launches a simulation with parameter multipliers set
        by p_vector and returns a single error value """
    
    # Open the analysis file
    with open(json_analysis_file_string, 'r') as f:
        orig_setup = json.load(f)
        
    # Copy the dict
    new_setup = copy.deepcopy(orig_setup)
    
    # Pull out the model_struct
    model_struct = orig_setup['FiberSim_setup']['model']
    
    # Deduce the base directory
    if (model_struct['relative_to'] == 'this_file'):
        model_base_dir = Path(json_analysis_file_string).parent.absolute()
    else:
        model_base_dir = model_struct['relative_to']

    # Set the working directory
    working_dir = os.path.join(model_base_dir,
                               model_struct['fitting']['working_folder'])
    
    # Clean the working directory
    working_dir = os.path.join(model_base_dir,
                               model_struct['fitting']['working_folder'])
    try:
        print('Trying to clean: %s' % working_dir)
        shutil.rmtree(working_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s' % (working_dir, e.strerror))
        
    # Check the working dir is there
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
        
    # Rename the fitting key
    new_setup['FiberSim_setup']['model']['manipulations'] = \
        new_setup['FiberSim_setup']['model'].pop('fitting')
        
    # Remove keys we don't need
    for k in ['working_folder', 'progress_folder','Python_objective_call']:
        new_setup['FiberSim_setup']['model']['manipulations'].pop(k)
        
    
        
    # Add Python objective call to characterization
    new_setup['FiberSim_setup']['characterization'][0] \
                ['post_sim_Python_call'] = \
        orig_setup['FiberSim_setup']['model']['fitting'] \
                ['Python_objective_call']
                
    # Write the new setup to file
    file_name = json_analysis_file_string.split('/')[-1]
    
    new_setup_file_string = os.path.join(working_dir,
                                         file_name)
    
    with open(new_setup_file_string, 'w') as f:
        json.dump(new_setup, f, indent=4)
    
    print(new_setup_file_string)
    print(p_vector)
    
    
    
        
    
        