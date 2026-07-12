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

        (char_setup_file, df_sample) = write_sample_setup_file(json_analysis_file_string,
                                                  sample_ind,
                                                  sample_values[sample_ind],
                                                  samples_folder,
                                                  chars_folder)

        # Save the char file
        char_files.append(char_setup_file)

        # Save the samples
        if (sample_ind == 0):
            all_samples = df_sample
        else:
            all_samples = pd.concat([all_samples, df_sample], ignore_index=True)

    # Save the sample data
    sample_data_file_string = str(Path(os.path.join(samples_folder, 'sample_values.csv')).resolve().absolute())
    print('Writing sample data to: %s' % sample_data_file_string)
    all_samples.to_csv(sample_data_file_string, index=False)

    return(char_files)
               
                
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

    # Hold the variable and multiplier as we go
    samples_dict = dict()
    
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

        # Now add the variable and the multiplier to the sampler dict
        if not ('kinetics' in adj['variable']):
            samples_dict[adj['variable']] = char_m
        else:
            var_name = ('%s_%i_%i_%i_%i' % (adj['variable'], adj['isotype'],adj['state'],adj['transition'],adj['parameter_number']))
            samples_dict[var_name] = char_m

    # Replace the adjustments with the new ones
    setup_dict['FiberSim_setup']['model']['manipulations']['adjustments'] = new_adj

    # Now we want to check for delta_hsl adjustments
    if ('delta_hsl' in setup_dict['FiberSim_setup']['model']['manipulations']):

        dhsl_dict = setup_dict['FiberSim_setup']['model']['manipulations']['delta_hsl']

        # We are going to need some samples
        no_of_dhsl_variables = len(dhsl_dict['data'])
        sampler = qmc.LatinHypercube(no_of_dhsl_variables, seed=sample_ind + 1)
        dhsl_values = sampler.random(1)
        
        # Loop through the characterizations
        for (char_ind, ch) in enumerate(setup_dict['FiberSim_setup']['characterization']):

            if (ch['type'] == 'twitch'):
                
            # delta_hsl stuff has to go in the protocol
                new_dhsl = dict()
                new_dhsl['type'] = dhsl_dict['type']

                for var_i in range(no_of_dhsl_variables):
                    span = dhsl_dict['data'][var_i]['factor_bounds'][1] - \
                            dhsl_dict['data'][var_i]['factor_bounds'][0]
                    sample_m = dhsl_values[0,var_i]
                    char_m = dhsl_dict['data'][var_i]['factor_bounds'][0] + (sample_m * span)
                    if ('factor_mode' in dhsl_dict['data'][var_i]) and \
                            (dhsl_dict['data'][var_i]['factor_mode'] == 'log'):
                        char_m = np.power(10, char_m)
                    new_dhsl[dhsl_dict['data'][var_i]['label']] = \
                        char_m * dhsl_dict['data'][var_i]['base_value']

                # Now set it
                # Get the protocol dict
                prot_dict = setup_dict['FiberSim_setup']['characterization'][char_ind]['protocol']

                # Add in the dhsl_dict
                for (pr_i, pr) in enumerate(prot_dict['data']):

                    pr['delta_hsl'] = new_dhsl
                    setup_dict['FiberSim_setup']['characterization'][char_ind]['protocol']['data'][pr_i] = pr

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

    if ('delta_hsl' in setup_dict['FiberSim_setup']['model']['manipulations']):
        del setup_dict['FiberSim_setup']['model']['manipulations']['delta_hsl']
    
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

    # Turn the samples_dict into a dataframe
    samples_dict
    df_samples = pd.DataFrame(data = samples_dict, index=[0])

    # Return the setup_file_string
    return (setup_file_string, df_samples) 


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
    