# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 12:20:56 2022

@author: kscamp3
"""

import os
import sys
import json
import copy

import numpy as np

from pathlib import Path

def characterize_models():
    """ Characterize models """
    
    # Detachment factor
    filament_compliance_factor = [1.0, 0.8, 0.6]
    
    # Variables
    base_setup_file_string = '../base/base_setup.json'
    base_model_file_string = '../base/base_model.json'
    
    generated_folder = '../generated'
    generated_setup_file_string = 'generated_setup.json'
    
    FiberSim_code_dir = '../../../../code/fiberpy/fiberpy'
    
    # Load the base setup and base model
    with open(base_setup_file_string, 'r') as f:
        base_setup = json.load(f) 

    with open(base_model_file_string, 'r') as f:
        base_model = json.load(f)
        
    # Get the parent dir, because all paths have to be converted to absolutes
    base_dir = Path(base_setup_file_string).parent.absolute()
        
    # Loop through mod factors
    rep_model_file_strings = list()
    
    for i in range(len(filament_compliance_factor)):

        # Copy the base model
        rep_model = copy.deepcopy(base_model)

        # Adjust the filament compliance
        # First the thin
        y = rep_model['thin_parameters']['a_k_stiff']
        rep_model['thin_parameters']['a_k_stiff'] = \
            y * filament_compliance_factor[i]
        
        # Now the thick
        y = rep_model['thick_parameters']['m_k_stiff']
        rep_model['thick_parameters']['m_k_stiff'] = \
            y * filament_compliance_factor[i]

        # Generate a model file name
        rfs = os.path.join(generated_folder,
                           'models',
                           'models_%i.json' % (i+1))
        
        # Correct the path
        rfs = str(Path(os.path.join(base_dir, rfs)).resolve())
        
        # Make sure the path exists
        parent_dir = Path(rfs).parent
        if not os.path.isdir(parent_dir):
            os.makedirs(parent_dir)

        # Write the model
        with open(rfs, 'w') as f:
            json.dump(rep_model, f, indent=4)

        # Add in to array
        rep_model_file_strings.append(rfs)
    
    # Now copy the setup
    generated_setup = copy.deepcopy(base_setup)
    
    # And change the file names to absolute paths
    generated_setup['FiberSim_characterization']['model']['relative_to'] = 'false'
    generated_setup['FiberSim_characterization']['model']['model_files'] = rep_model_file_strings
    generated_setup['FiberSim_characterization']['model']['options_file'] = \
        str(Path(os.path.join(base_dir,
                              generated_setup['FiberSim_characterization']['model']['options_file'])).resolve())
        
    # Load up the characterizations 
    characterize = generated_setup['FiberSim_characterization']['characterization']
    
    # Loop through the characterizations, changing these paths
    for i,c in enumerate(characterize):
        c['relative_to'] = 'this_file'
        c['sim_folder'] = str(Path(os.path.join(base_dir,
                                                c['sim_folder'])).resolve())
        generated_setup['FiberSim_characterization']['characterization'][i] = c

    # And finally the setup file
    generated_setup_file_string = str(Path(os.path.join(base_dir,
                                                        generated_folder,
                                                        generated_setup_file_string)).resolve())
    with open(generated_setup_file_string, 'w') as f:
        json.dump(generated_setup, f, indent=4)
        
    # Generate a command line
    cs = 'pushd \"%s\" & python FiberPy.py characterize %s & popd' % \
            (FiberSim_code_dir, generated_setup_file_string)
    
    # And run it
    os.system(cs)
    

if __name__ == "__main__":
    characterize_models()
