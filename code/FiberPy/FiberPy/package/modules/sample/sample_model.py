# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 15:35:50 2023

@author: Campbell
"""

import os
import json
import shutil
import copy

import multiprocessing
import threading
import subprocess

import numpy as np
import pandas as pd

from pathlib import Path

from scipy.stats import qmc

from ..protocols import protocols as prot
from ..batch import batch

from ..characterize.characterize_utilities import \
    return_FiberCpp_exe_dict, \
    return_base_dir, \
    return_model_file_strings, \
    return_options_file_string, \
    prepare_clean_dir

def sample_model(json_analysis_file_string):
    """ Code takes a json struct that includes a model file
        and samples the model """

    # Check for the analysis file
    if not (json_analysis_file_string):
        print('sample_model: no analysis file specified. Exiting')
        exit(1)
        
    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        anal_struct = json_data['FiberSim_setup']
      
    # Check that we have a sampling model
    if not ("sampling" in anal_struct["model"]):
        print('Error, no sampling structure in: %s' %
              json_analysis_file_string)
        exit(1)

    # Generate the characterization files
    characterization_files = \
        generate_characterization_files(json_analysis_file_string)

    # Generate a sequence of command strings
    command_strings = []
        
    for cf in characterization_files:
        print(cf)
        cs = 'python FiberPy.py characterize %s' % cf
        command_strings.append(cs)

    # Now run them        
    batch_command_strings(command_strings)

    # Run post-Python_function
    if ('sampling' in anal_struct['model']):
        if ('post_sim_Python_call' in anal_struct['model']['sampling']):
            python_file = anal_struct['model']['sampling']['post_sim_Python_call']
            
            if (anal_struct['model']['relative_to'] == 'this_file'):
                working_dir = Path(json_analysis_file_string).parent.absolute()
            else:
                working_dir = ''
                
            command_string = 'python %s' % \
                os.path.join(working_dir, python_file)

            subprocess.call(command_string)   


def generate_characterization_files(json_analysis_file_string):
    """ Generates a sequence of characterization files to sample a model
        over a defined parameter space """
    
    # First load the sampling file and pull off the relevant structs
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        model_struct = json_data['FiberSim_setup']['model']
        sampling_struct = model_struct['sampling']
        char_struct = json_data['FiberSim_setup']['characterization']
    
    # Deduce the base model file string
    base_dir = return_base_dir(json_analysis_file_string, 'model')
    base_model_file_string = str(Path(os.path.join(base_dir,
                                                   sampling_struct['base_model'])).resolve().absolute())
    # Deduce the base options file string
    base_options_file_string = return_options_file_string(json_analysis_file_string)

    # We are going to generate a sequence of setup files, one for each model sample
    # They have to go into a clean directory that persists until the simulations are complete
    # So we will make a top folder, clean that, then put the setup_files in a sub-folder called 'samples'
    # that uses a parallel-sub-folder called 'chars' for the working space

    # top_folder
    top_folder = str(Path(os.path.join(base_dir,
                                       sampling_struct['generated_folder'])).resolve().absolute())
    prepare_clean_dir(top_folder)

    # samples_folder
    samples_folder = str(Path(os.path.join(top_folder, 'samples')).resolve().absolute())
    prepare_clean_dir(samples_folder)

    # chars_folder
    chars_folder = str(Path(os.path.join(top_folder, 'chars')).resolve().absolute())
    prepare_clean_dir(chars_folder)

    # Now deduce parameters for the sampling
    adjustments = sampling_struct['adjustments']
    no_of_parameters = len(adjustments)
    no_of_samples = sampling_struct['no_of_samples']

    # Generate the sample values
    sampler = qmc.LatinHypercube(no_of_parameters, seed=1)
    sample_values = sampler.random(no_of_samples)
    
    print('no_of_parameters: %i' % no_of_parameters)
    print('no_of_samples: %i' % no_of_samples)
    print(sample_values)

    # Loop through the samples, storing them as we go

    char_files = []

    for sample_ind in range(no_of_samples):

        char_setup_file = write_sample_setup_file(json_analysis_file_string,
                                                  sample_ind,
                                                  sample_values[sample_ind],
                                                  samples_folder,
                                                  chars_folder)

        char_files.append(char_setup_file)

    return(char_files)

        
    # # We also need to wipe the char folder
    # generated_char_dir = '%s_char' % top_generated_dir
    # prepare_clean_dir(generated_char_dir)
        
    # # Finally, we need to prep the sim directory for each charactization
    # # Keep track of the characterization directories as we go
    # ch_dirs = [];
    # sim_dirs = []
    # for (char_id, ch) in enumerate(char_struct):

    #     char_dir = return_base_dir(json_analysis_file_string,
    #                                'characterization',
    #                                dict_index = char_id)

    #     ch_dirs.append(char_dir)

    #     sim_dir_id = str(Path(os.path.join(char_dir,
    #                                        ch['sim_folder'])).resolve().absolute())
    #     prepare_clean_dir(sim_dir_id)

    #     sim_dirs.append(sim_dir_id)
        
    # # Now deduce parameters for the sampling
    # adjustments = sampling_struct['adjustments']
    # no_of_parameters = len(adjustments)
    # no_of_samples = sampling_struct['no_of_samples']

    # # Generate the sample values
    # sampler = qmc.LatinHypercube(no_of_parameters, seed=1)
    # sample_values = sampler.random(no_of_samples)
    
    # print('no_of_parameters: %i' % no_of_parameters)
    # print('no_of_samples: %i' % no_of_samples)
    # print(sample_values)

    # # Load up the base model
    # with open(base_model_file_string, 'r') as f:
    #     base_model = json.load(f)
        
    # # Prepare a list of characterization files
    # characterization_file_strings = []
       
    # for sample_counter in range(no_of_samples):
        
    #     # Copy the setup_up struct and convert to a characterization file
    #     sample_characterize = copy.deepcopy(json_data)
        
    #     # Create the gen directory for the sample
    #     sample_gen_dir = os.path.join(top_generated_dir,
    #                                   ('sample_%i' % (sample_counter + 1)))
        
    #     prepare_clean_dir(sample_gen_dir)

    #     # Update the FiberCpp_exe section
    #     sample_characterize['FiberSim_setup']['FiberCpp_exe'] = \
    #         return_FiberCpp_exe_dict(json_analysis_file_string)
             
    #     # Update the model section
    #     sample_characterize['FiberSim_setup']['model']['relative_to'] = \
    #         'False';
        
    #     # Set and copy the options file
    #     temp, options_file_end = os.path.split(base_options_file_string)
    #     new_options_file_string = str(Path(
    #         os.path.join(sample_gen_dir, options_file_end)).absolute().resolve())
    #     sample_characterize['FiberSim_setup']['model']['options_file'] = \
    #         new_options_file_string
            
    #     shutil.copy(base_options_file_string,
    #                 new_options_file_string)
        
    #     # Delete the sampling and replace with manipulations
    #     del sample_characterize['FiberSim_setup']['model']['sampling']

    #     # Create the manipulations
    #     sample_characterize['FiberSim_setup']['model'] \
    #         ['manipulations'] = dict()
        
    #     # Create and copy the model file
    #     temp, model_file_end = os.path.split(base_model_file_string)
    #     new_model_file_string = str(Path(os.path.join(sample_gen_dir,
    #                                          model_file_end)).absolute().resolve())
    #     sample_characterize['FiberSim_setup']['model'] \
    #         ['manipulations']['base_model'] = new_model_file_string
            
    #     shutil.copy(base_model_file_string, new_model_file_string)
        
    #     # Set the generated fold
    #     # This path has to be adjusted a bit so that it is different from the
    #     # sampling folder, and thus won't be wiped at the beginning of the
    #     # simulations
    #     temp, last_bit = os.path.split(sample_gen_dir)
    #     temp_generated_dir = os.path.join(('%s_char' % temp),
    #                                       last_bit)
        
    #     sample_characterize['FiberSim_setup']['model'] \
    #         ['manipulations']['generated_folder'] = temp_generated_dir

    #     # Now loop through the characterize structs - we need to do this here
    #     # in case we have to adapt the twitch protocol for the characterization
    #     for (char_id, ch) in enumerate(sample_characterize['FiberSim_setup']['characterization']):
            
    #         if (ch['type'] == 'twitch'):
    #             if ('data' in ch['protocol']):
    #                 tw_prot_data = ch['protocol']['data'][0]

    #         # Make an array of adjustments
    #         adjusts = []
        
    #         # Make a dict of parameter values
    #         par_values = dict()
        
    #         for par_counter in range(no_of_parameters):
    #             # Copy the adjustments
    #             sample_adj = sampling_struct['adjustments'][par_counter]
    #             characterize_adj = copy.deepcopy(sample_adj)
            
    #             if (sample_adj['variable'].startswith('Ca_transient_')):
    #                 # Pull off the key
    #                 twitch_key = sample_adj['variable'].split('Ca_transient_')[-1]
                
    #                 # Get the base value
    #                 base_value = tw_prot_data[twitch_key]
                
    #                 # Now deduce the multiplier
    #                 span = sample_adj['factor_bounds'][1] - sample_adj['factor_bounds'][0]
                
    #                 sample_m = sample_values[sample_counter][par_counter]
                
    #                 characterize_m = sample_adj['factor_bounds'][0] + \
    #                     (sample_m * span)
                    
    #                 if ('factor_mode' in sample_adj) and (sample_adj['factor_mode'] == 'log'):
    #                     characterize_m = np.power(10, characterize_m)
                    
    #                 tw_prot_data[twitch_key] = characterize_m * base_value
                
    #                 # Store the value
    #                 par_values[sample_adj['variable']] = tw_prot_data[twitch_key]
                
    #                 continue
            
    #             if ((sample_adj['variable'] == 'm_kinetics') or
    #                     (sample_adj['variable'] == 'c_kinetics')):

    #                 # Special case for kinetics
    #                 if ('extension' in sample_adj):
    #                     base_value = base_model[sample_adj['variable']][sample_adj['isotype']-1]\
    #                         ['state'][sample_adj['state']-1]['extension']
                        
    #                     # Store the key
    #                     par_key = '%s_isotype_%i_state_%i_extension' % \
    #                                 (sample_adj['variable'], sample_adj['isotype'],
    #                                  sample_adj['state'])
    #                 else:
    #                     # Transition parameters
    #                     y = np.asarray(base_model[sample_adj['variable']][sample_adj['isotype']-1] \
    #                                        ['state'][sample_adj['state']-1] \
    #                                        ['transition'][sample_adj['transition']-1]\
    #                                        ['rate_parameters'],
    #                               dtype = np.float32)
                    
    #                     base_value = y[sample_adj['parameter_number'] - 1]
                    
    #                     # Store the key
    #                     par_key = '%s_isotype_%i_scheme_%i_transition_%i_parameter_%i' % \
    #                                 (sample_adj['variable'], sample_adj['isotype'],
    #                                  sample_adj['state'], sample_adj['transition'],
    #                                  sample_adj['parameter_number'])
                
    #                 # Now work out the values
    #                 characterize_adj['base_value'] = float(base_value)
                
    #                 # Now deduce the multiplier
    #                 span = sample_adj['factor_bounds'][1] - sample_adj['factor_bounds'][0]
                
    #                 sample_m = sample_values[sample_counter][par_counter]
                
    #                 characterize_m = sample_adj['factor_bounds'][0] + \
    #                     (sample_m * span)
                    
    #                 if ('factor_mode' in sample_adj) and (sample_adj['factor_mode'] == 'log'):
    #                     characterize_m = np.power(10, characterize_m)

    #                 characterize_adj['multipliers'] = []
    #                 characterize_adj['multipliers'].append(characterize_m)
                
    #                 characterize_adj['output_type'] = 'float'
                
    #                 # Add it in
    #                 adjusts.append(characterize_adj)
                               
    #                 par_values[par_key] = base_value * characterize_m
                
    #                 continue
            
            
    #             # Everything else
    #             base_value = base_model[sample_adj['class']][sample_adj['variable']]
    #             span = sample_adj['factor_bounds'][1] - sample_adj['factor_bounds'][0]
    #             sample_m = sample_values[sample_counter][par_counter]
    #             characterize_m = sample_adj['factor_bounds'][0] + \
    #                 (sample_m * span)
    #             if ('factor_mode' in sample_adj) and (sample_adj['factor_mode'] == 'log'):
    #                 characterize_m = np.power(10, characterize_m)
    #             characterize_adj['multipliers'] = []
    #             characterize_adj['multipliers'].append(characterize_m)
    #             characterize_adj['output_type'] = 'float'

    #              # Add it in
    #             adjusts.append(characterize_adj)
            
    #             # Store the value
    #             par_key = '%s_%s' % (sample_adj['class'], sample_adj['variable'])
    #             par_values[par_key] = base_value * characterize_m
                
    #         # Make a dataframe from the par_values
    #         par_df = pd.DataFrame([par_values])
        
    #         if (sample_counter == 0):
    #             collated_parameters = par_df
    #         else:
    #             collated_parameters = pd.concat([collated_parameters, par_df],
    #                                             ignore_index = True)
        
    #         # Add the adjustments into sample_characterize
    #         sample_characterize['FiberSim_setup'] \
    #             ['model']['manipulations']['adjustments'] = adjusts
                
    #         if (ch['type'] == 'unloaded_shortening'):
    #             new_ch = characterize_unloaded_shortening(json_analysis_file_string,
    #                                                   sample_counter,
    #                                                   sample_gen_dir,
    #                                                   tw_protocol)
        
    #         else:
    #             # Adjust the output dir
    #             new_ch = copy.deepcopy(ch)

    #             new_ch['relative_to'] = 'False'
    #             new_ch['sim_folder'] = str(Path(os.path.join(sim_dirs[char_id],
    #                                                          ('sample_%i' % (sample_counter + 1)))).
    #                                         absolute().resolve())

    #             if ('protocol' in ch):
    #                 if ('protocol_folder' in ch['protocol']):
    #                     new_ch['protocol']['protocol_folder'] = \
    #                         str(Path(os.path.join(new_ch['sim_folder'], ch['protocol']['protocol_folder'])).
    #                             absolute().resolve())
           
    #         # Adjust the post-sim Python call
    #         if ('post_sim_Python_call' in ch):
    #             new_ch['post_sim_Python_call'] = str(Path(
    #                 os.path.join(char_dirs(ch_id),ch['post_sim_Python_call'])).resolve().absolute())

    #         # Repack
    #         sample_characterize['FiberSim_setup']['characterization'][char_id] = new_ch
        
    #         # Create a file name for the characterization file
    #         characterization_file_string = \
    #             os.path.join(sample_gen_dir,
    #                          ('characterize_%i.json' % (sample_counter+1)))
            
    #         with open(characterization_file_string, 'w') as f:
    #             json.dump(sample_characterize, f, indent=4)
            
    #         # Append to the list
    #         characterization_file_strings.append(characterization_file_string)
        
    # # Display the parameter values
    # print(collated_parameters)
    
    # # Output to file
    # par_file_string = os.path.join(top_generated_dir, 'parameter_values.xlsx')
    # collated_parameters.to_excel(par_file_string, index=False)
        
    # # Return
    # return characterization_file_strings        
                
                
def characterize_unloaded_shortening(json_analysis_file_string,
                                     sample_counter,
                                     sample_generated_directory,
                                     tw_protocol):
    """ Adjusts the characterization dict to handle unloaded shortening """
    
    # First load the file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        
    # Replace the twitch protocol
    json_data['FiberSim_setup']['characterization'][0]['twitch_protocol'] = tw_protocol
    
    # Pull off the characterization component
    orig_struct = json_data['FiberSim_setup']['characterization'][0]
    
    # Copy it
    new_struct = copy.deepcopy(orig_struct)
    
    new_struct['type'] = 'freeform'
    new_struct['sim_folder'] = str(Path(os.path.join('..',
                                            orig_struct['sim_folder'],
                                            ('sample_%i' % (sample_counter + 1)))).resolve())
    
    # Now deduce some folders
    if (orig_struct['relative_to'] == 'this_file'):
        parent_dir = Path(json_analysis_file_string).parent.absolute()
        
        new_struct['relative_to'] = 'False'
        new_struct['sim_folder'] = os.path.join(str(
                Path(parent_dir,
                     orig_struct['sim_folder'],
                     ('sample_%i' % (sample_counter + 1)))))
        
        # Make a protocol file string
        protocol_file_string = os.path.join(sample_generated_directory,
                                        ('protocol_%i.txt' % (sample_counter + 1)))
    else:
        print('More work on paths required')
        exit(1)
        
    # Now make the protocol
    pr_df = prot.create_twitch_protocol(
                time_step = new_struct['twitch_protocol']['time_step_s'],
                n_points = new_struct['twitch_protocol']['n_points'],
                stimulus_times_s = new_struct['twitch_protocol']['stimulus_time_s'],
                Ca_content = new_struct['twitch_protocol']['Ca_content'],
                stimulus_duration_s = new_struct['twitch_protocol']['stimulus_duration_s'],
                k_leak = new_struct['twitch_protocol']['k_leak'],
                k_act = new_struct['twitch_protocol']['k_act'],
                k_serca = new_struct['twitch_protocol']['k_serca'],
                mode_vector = np.zeros(new_struct['twitch_protocol']['n_points']))
    
    # Write the protocol to file
    prot.write_protocol_to_file(pr_df, protocol_file_string)
    
    # Add the protocol_file_string to the char struct
    new_struct['protocol_files'] = []
    new_struct['protocol_files'].append(protocol_file_string)
    
    # Delete the twitch data
    del new_struct['twitch_protocol']
    
    # Return the struct
    return new_struct

        # sample_setup_file = write_sample_setup_file(json_analysis_file,
        #                                             sample_counter,
        #                                             sample_values[sample_counter],
        #                                             samples_folder)

def write_sample_setup_file(json_analysis_file_string,
                            sample_ind,
                            sample_values,
                            samples_folder,
                            chars_folder):
    """ Writes a sample setup file for a given sample """

    # Loads the original setup file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    # Copy it
    setup_dict = copy.deepcopy(json_dict)

    # Make sure paths are fixed
    setup_dict['FiberSim_setup']['FiberCpp_exe'] = return_FiberCpp_exe_dict(json_analysis_file_string)

    # Adjust the options file
    setup_dict['FiberSim_setup']['model']['relative_to'] = 'False'
    setup_dict['FiberSim_setup']['model']['options_file'] = \
        str(Path(return_options_file_string(json_analysis_file_string)).resolve().absolute())

    # Copy the sampling_dict to an maipulations_dict
    setup_dict['FiberSim_setup']['model']['manipulations'] = \
        copy.deepcopy(setup_dict['FiberSim_setup']['model']['sampling'])

    # Adjust the base_model
    base_dir = return_base_dir(json_analysis_file_string,
                               dict_key = 'model')
    base_model_file_string = str(Path(os.path.join(base_dir,
                                                   setup_dict['FiberSim_setup']['model']['manipulations']['base_model'])).
                                 resolve().absolute())
    setup_dict['FiberSim_setup']['model']['manipulations']['base_model'] = base_model_file_string

    # Adjust the generated folder
    generated_folder = os.path.join(chars_folder, ('%i' % (sample_ind + 1)))
    setup_dict['FiberSim_setup']['model']['manipulations']['generated_folder'] = \
        generated_folder

    # Now span through the adjustments, creating an array of new ones
    new_adj = []
    
    for (adj_ind, adj) in enumerate(setup_dict['FiberSim_setup']['model']['manipulations']['adjustments']):

        # If the manipulation is not to do with a Ca transient, we can store it as a new adjustment
        # and handle the effect in the setup file
        if (not adj['variable'].startswith('Ca_transient')):

            span = adj['factor_bounds'][1] - adj['factor_bounds'][0]
            sample_m = sample_values[adj_ind]
            char_m = adj['factor_bounds'][0] + (sample_m * span)

            if ('factor_mode' in adj) and (adj['factor_mode'] == 'log'):
                char_m = np.power(10, char_m)

            # Append it
            a = copy.deepcopy(adj)
            a['multipliers'] = [char_m]
            a['output_type'] = 'float'

            new_adj.append(a)

        else:
            # We have to go into the protocol and adjust the Ca transient values directly
            for (char_ind, ch) in enumerate(setup_dict['FiberSim_setup']['characterization']):

                if (not ('protocol' in ch)):
                    print('Error: trying to adjust a Ca_transient_ variable in characterization[%i] but no protocol found. Exiting' % char_ind)
                    exit(1)

                if (not ('data' in ch['protocol'])):
                    print('Error: trying to adjust a Ca_transient_ variable in characterization[%i] but no data found. Exiting' % char_ind)
                    exit(1)

                for (data_ind, d) in enumerate(ch['protocol']['data']):

                    twitch_key = adj['variable'].split('Ca_transient_')[-1]
                    base_value = ch['protocol']['data'][0][twitch_key]
                    span = adj['factor_bounds'][1] - adj['factor_bounds'][0]
                    sample_m = sample_values[adj_ind]
                    char_m = adj['factor_bounds'][0] + (sample_m * span)
                    if ('factor_mode' in adj) and (adj['factor_mode'] == 'log'):
                        char_m = np.power(10, char_m)
                    setup_dict['FiberSim_setup']['characterization'][char_ind]['protocol']['data'][0][twitch_key] = \
                        char_m * base_value

    # Repalce the adjustments with the new ones
    setup_dict['FiberSim_setup']['model']['manipulations']['adjustments'] = new_adj

    # Getting there, clean up the folders for the characterization
    for (char_ind, ch) in enumerate(setup_dict['FiberSim_setup']['characterization']):
        
        setup_dict['FiberSim_setup']['characterization'][char_ind]['relative_to'] = 'False'
        
        base_dir = return_base_dir(json_analysis_file_string,
                                   'characterization',
                                   dict_index = char_ind)

        sim_folder = str(Path(os.path.join(base_dir,
                                           ch['sim_folder'],
                                           ('%i' % (sample_ind + 1)))).resolve().absolute())
        setup_dict['FiberSim_setup']['characterization'][char_ind]['sim_folder'] = sim_folder

        if ('protocol' in ch) and ('protocol_folder' in ch['protocol']):
            prot_folder = str(Path(os.path.join(sim_folder, ch['protocol']['protocol_folder'])).resolve().absolute())
            setup_dict['FiberSim_setup']['characterization'][char_ind]['protocol']['protocol_folder'] = prot_folder

    # Delete some stuff we don't need
    del setup_dict['FiberSim_setup']['model']['manipulations']['no_of_samples']
    del setup_dict['FiberSim_setup']['model']['sampling']
    
    # Create the folder for the setup file
    sample_folder = os.path.join(samples_folder, ('%i' % (sample_ind + 1)))
    prepare_clean_dir(sample_folder)

    # Deduce the name
    setup_file_string = str(Path(os.path.join(sample_folder,
                                              str(Path(json_analysis_file_string).name))).
                            absolute().resolve())

    # Write it
    print('Writing sample setup file: %s' % setup_file_string)
    with open(setup_file_string, 'w') as f:
        json.dump(setup_dict, f, indent=4)

    # Return the setup_file_string
    return setup_file_string


def batch_command_strings(command_strings, figures_only=False):
    """ Runs a list of command strings as separate threads """
    
    if not figures_only:
        my_list = command_strings
            
        # Get max threads available
        num_processes = multiprocessing.cpu_count()-1
        print('Running batch using %i threads' % num_processes)
        
        threads = []
        while threads or my_list:
            if (len(threads) < num_processes) and my_list:
                t = threading.Thread(target=worker, args=[my_list.pop()])
                t.setDaemon(True)
                t.start()
                threads.append(t)
            else:
                for thread in threads:
                    if not thread.is_alive():
                        threads.remove(thread)
                        
def worker(cmd):
    subprocess.call(cmd)
    