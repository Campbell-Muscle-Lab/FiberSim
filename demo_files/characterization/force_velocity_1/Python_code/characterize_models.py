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
    detachment_factor = [1.0, 0.5, 2.0]
    
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
    
    for i in range(len(detachment_factor)):

        # Copy the base model
        rep_model = copy.deepcopy(base_model)

        # Detachment rate
        # Pull off the transition vector
        y = np.asarray(rep_model['m_kinetics'][0]['scheme'][2]['transition'][0]['rate_parameters'],
                       dtype=np.float32)
        y[0] = y[0] * detachment_factor[i]
         # Rewrite array into json format
        rep_model['m_kinetics'][0]['scheme'][2]['transition'][0]['rate_parameters'] = \
            y.tolist()

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
    
    # # Find the freeform one and add in protocol files
    # for i, c in enumerate(characterize):
    #     if not (c['type'] == 'freeform'):
    #         continue
        
    #     prot_points = round(c['sim_duration_s'] / c['sim_time_step_s'])
    #     delta_hsl = np.zeros(prot_points)
    #     t = np.cumsum(c['sim_time_step_s'] * np.ones(prot_points))
    #     step_ind = np.argmax(t > c['freeform_step_time_s'])
        
    #     step_size_nm = c['freeform_rel_step_size'] * base_model['muscle']['initial_hs_length']
    #     delta_hsl[step_ind] = step_size_nm
        
    #     prot_file_strings = []
        
    #     for j,p in enumerate(c['freeform_pCa_values']):
    #         df = prot.create_length_control_protocol(
    #                 time_step=c['sim_time_step_s'],
    #                 n_points=prot_points,
    #                 initial_pCa=8.0,
    #                 step_pCa=p,
    #                 delta_hsl=delta_hsl)
            
    #         # Write to file
    #         prot_folder = os.path.join(base_dir,
    #                                    generated_folder,
    #                                    'protocols')
            
    #         if not os.path.isdir(prot_folder):
    #             os.makedirs(prot_folder)
            
    #         prot_file_string = os.path.join(prot_folder,
    #                                         ('%s_%i.txt' % (c['prot_file_base'], j)))
    #         prot_file_string = str(Path(prot_file_string).absolute().resolve())
    #         prot_file_strings.append(prot_file_string)
    #         df.to_csv(Path(prot_file_string).absolute().resolve(),
    #                   sep='\t',
    #                   index=None)
            
    #     characterize[i]['protocol_files'] = prot_file_strings
    
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
