
import os
import json

import shutil

import numpy as np
import pandas as pd

from pathlib import Path

from .characterize_utilities import \
        create_sim_input_and_output_dirs, \
        prepare_clean_dir, \
        prepare_protocols, \
        prepare_repeats, \
        prepare_simulation_dir, \
        return_base_dir, \
        return_batch_figs_dict, \
        return_hs_lengths, \
        return_FiberCpp_exe_dict, \
        return_model_file_strings, \
        return_options_file_string, \
        return_run_mode, \
        update_and_write_model_file, \
        update_setup_file_string_with_new_char

from ..batch import batch

from ..utilities import utilities as ut

def twitch_loads(json_analysis_file_string, char_index = 0):
    """ Twitch loads runs simulations of twitch contractions where the muscle shortens
    against defined loads that can be imposed:
    1) at a fixed time
    2) at a fixed time offset from peak force
    3) as soon as isometric force reaches the load
    """
    run_mode = return_run_mode(json_analysis_file_string,
                               char_index = char_index)

    model_file_strings = return_model_file_strings(json_analysis_file_string)

    hs_lengths = return_hs_lengths(json_analysis_file_string,
                                   char_index = char_index)

    top_sim_dir = return_base_dir(json_analysis_file_string,
                                 'characterization',
                                 append_key = 'sim_folder',
                                 dict_index = char_index)

    # Prepare the top sim dir
    prepare_simulation_dir(top_sim_dir, run_mode)

    # Pull off the char_dict
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)
        char_dict = json_dict['FiberSim_setup']['characterization'][char_index]

    # Check for protocols, which can be used for twitches
    # If not present, the protocol loop will just run once
    if ('protocol' in char_dict):
        protocol_file_strings = prepare_protocols(json_analysis_file_string,
                                                  char_index = char_index)

        no_of_protocols = len(protocol_file_strings)
    else:
        protocol_file_strings = []
        no_of_protocols = 1

    # Create a working version of the setup file - this will be modified
    # in successive loops to implement different modes
    working_dir = return_base_dir(json_analysis_file_string,
                        'characterization',
                        append_key = 'working_folder',
                        dict_index = char_index)
    
    prepare_clean_dir(working_dir)
    
    working_analysis_file_string = os.path.join(
        working_dir,
        str(Path(json_analysis_file_string).name) )

    # Adjust the model paths to account for the new working file
    setup_folder = str(Path(json_analysis_file_string).parent.name)
    json_dict['FiberSim_setup']['model']['options_file'] = \
        os.path.join('..',
                     setup_folder,
                     json_dict['FiberSim_setup']['model']['options_file'])
    for (i, mf) in enumerate(json_dict['FiberSim_setup']['model']['model_files']):
        json_dict['FiberSim_setup']['model']['model_files'][i] = \
            os.path.join('..',
                         setup_folder,
                         mf)
    
    # Tidy up
    working_analysis_file_string = str(Path(working_analysis_file_string).resolve())

    # Set a default mode
    if not ('trial_modes' in char_dict):
        char_dict['trial_modes'] = ['isometric']

    # Set the number of modes
    no_of_trial_modes = len(char_dict['trial_modes'])

    # Prepare the directories
    mode_sim_dir = []
    for i in range(no_of_trial_modes):
        mode_sim_dir.append(os.path.join(top_sim_dir, ('%i' % (i+1))))
        prepare_simulation_dir(mode_sim_dir[i], run_mode)

    # Loop through the simulations twice, the first time
    # we are in isometric mode to establish the peak time and force
    # the second time, we are in loaded mode

    for mode_i in range(no_of_trial_modes):

        # Write the setup_file - first time around it's the original
        # analysis file, after that, it's being modified
        update_setup_file_string_with_new_char(
            json_dict, char_dict, char_index,
            working_analysis_file_string)

        # Set some stuff up
        sim_dir = mode_sim_dir[mode_i]
        if (char_dict['trial_modes'][mode_i] == 'isometric'):
            impose_afterloads = False
        else:
            impose_afterloads = True

        # Set up a dir counter
        dir_counter = 0;

        # Set up a jobs array to hold the repeat jobs
        all_jobs = []

        # Loop through model files and hs_lengths
        for (mod_i, mfs) in enumerate(model_file_strings):

            for (hsl_i, hsl) in enumerate(hs_lengths):

                # Update counter
                dir_counter = dir_counter + 1

                # Create sim_input and sim_output folders
                (sim_input_dir, sim_output_dir) = \
                    create_sim_input_and_output_dirs(sim_dir,
                                                     dir_counter = dir_counter)

                print(sim_output_dir)
                print(hs_lengths)

                # Set the new model_file_string
                new_model_file_string= os.path.join(sim_input_dir,
                                                   Path(mfs).name)
                # Tidy it
                new_model_file_string = str(Path(new_model_file_string).resolve())

                # Write it
                update_and_write_model_file(mfs,
                                            new_model_file_string,
                                            char_dict,
                                            hsl)

                # Loop through protocols
                for prot_index in range(no_of_protocols):

                    # Write the protocol if you need to
                    if (protocol_file_strings):
                        # We need to copy one
                        file_name = Path(protocol_file_strings[prot_index]).name
                        new_protocol_file_string = os.path.join(sim_input_dir, file_name)
                        shutil.copy(protocol_file_strings[prot_index], new_protocol_file_string)

                    # Create needed repeats
                    repeat_jobs = prepare_repeats(working_analysis_file_string,
                                                  sim_input_dir,
                                                  sim_output_dir,
                                                  new_model_file_string,
                                                  protocol_file_strings[prot_index],
                                                  char_index = char_index,
                                                  prot_index = prot_index,
                                                  impose_afterloads = impose_afterloads)

                    # Append the jobs
                    for j in repeat_jobs:
                        all_jobs.append(j)

        # Create a batch

        # We will need some folders
        sim_output_dir = os.path.join(sim_dir, 'sim_output')

        isometric_batch = dict()
        isometric_batch['FiberSim_batch'] = dict()
        isometric_batch['FiberSim_batch']['FiberCpp_exe'] = \
            return_FiberCpp_exe_dict(json_analysis_file_string)
        isometric_batch['FiberSim_batch']['job'] = all_jobs
        isometric_batch['FiberSim_batch']['batch_figures'] = \
            return_batch_figs_dict(json_analysis_file_string,
                                   sim_output_dir,
                                   char_index)

        # Write the batch
        batch_file_string = os.path.join(sim_dir, 'batch.json')

        with open(batch_file_string, 'w') as f:
            json.dump(isometric_batch, f, indent = 4)

        # Now run it
        if not (run_mode == 'figures_only'):
            batch.run_batch(batch_file_string)

        # If there is going to be another loop, analyze the first mode (assumed isometric)
        if ( (mode_i == 0) and (no_of_trial_modes > 1) ):
            twitch_data = analyze_twitches(mode_sim_dir[mode_i])

        # Now update the char_dict for different modes
        next_mode_i = mode_i + 1

        if (next_mode_i < no_of_trial_modes):

            if (char_dict['trial_modes'][next_mode_i] == 'afterload'):
                rel_loads = char_dict['protocol']['data'][prot_index]['afterload']['rel_load']
                for i in range(len(rel_loads)):
                    char_dict['protocol']['data'][prot_index]['afterload']['load'][i] = \
                        twitch_data['pas_force'] + \
                            (rel_loads[i] * (twitch_data['max_force'] - 
                                             twitch_data['pas_force']))
                    char_dict['protocol']['data'][prot_index]['afterload']['min_init_time_s'][i] = 0.0

            if (char_dict['trial_modes'][next_mode_i] == 'release_at_max_force'):
                rel_loads = char_dict['protocol']['data'][prot_index]['afterload']['rel_load']
                for i in range(len(rel_loads)):
                    char_dict['protocol']['data'][prot_index]['afterload']['load'][i] = \
                        twitch_data['pas_force'] + \
                            (rel_loads[i] * (twitch_data['max_force'] - 
                                             twitch_data['pas_force']))
                    char_dict['protocol']['data'][prot_index]['afterload']['min_init_time_s'][i] = \
                        twitch_data['time_s_at_max_force']

            if (char_dict['trial_modes'][next_mode_i] == 'release_at_max_thin_act'):
                rel_loads = char_dict['protocol']['data'][prot_index]['afterload']['rel_load']
                for i in range(len(rel_loads)):
                    char_dict['protocol']['data'][prot_index]['afterload']['load'][i] = \
                        twitch_data['pas_force'] + \
                            (rel_loads[i] * (twitch_data['max_force'] - 
                                             twitch_data['pas_force']))
                    char_dict['protocol']['data'][prot_index]['afterload']['min_init_time_s'][i] = \
                        twitch_data['time_s_at_max_thin_act']

            if (char_dict['trial_modes'][next_mode_i] == 'release_at_defined_time'):
                rel_loads = char_dict['protocol']['data'][prot_index]['afterload']['rel_load']
                for i in range(len(rel_loads)):
                    char_dict['protocol']['data'][prot_index]['afterload']['load'][i] = \
                        twitch_data['pas_force'] + \
                            (rel_loads[i] * (twitch_data['max_force'] - 
                                             twitch_data['pas_force']))

def analyze_twitches(sim_dir, twitch_pas_points = 90):
    """ Analyzes the twitches to get the peak force and times """

    # Pull off the results files
    results_files = ut.return_sim_results_files_in_nested_dir(sim_dir);

    # Prepare arrays for max force and time
    no_of_files = len(results_files)
    max_force_array = np.nan * np.ones(no_of_files)
    max_force_time_s_array = np.nan * np.ones(no_of_files)
    max_thin_act_array = np.nan * np.ones(no_of_files)
    max_thin_act_time_s_array = np.nan * np.ones(no_of_files)
    pas_force_array = np.nan * np.ones(no_of_files)
    pas_force_time_s_array = np.nan * np.ones(no_of_files)

    # Cycle through the files
    for (i, f) in enumerate(results_files):

        d = pd.read_csv(f, sep='\t')

        max_force_array[i] = d['m_force'].max()
        max_force_idx = d['m_force'].idxmax()
        max_force_time_s_array[i] = d['time'].iloc[max_force_idx]
        max_thin_act_array[i] = d['hs_1_a_pop_2'].max()
        max_thin_act_idx = d['hs_1_a_pop_2'].idxmax()
        max_thin_act_time_s_array[i] = d['time'].iloc[max_thin_act_idx]

        # Calculate the first derivative of force to try and get twitch start
        d['dm_force_dt'] = d['m_force'].diff()
        max_dfdt_idx = d['dm_force_dt'].idxmax()
        pas_force_idx = twitch_pas_points
        pas_force_array[i] = d['m_force'].iloc[twitch_pas_points]
        pas_force_time_s_array[i] = d['time'].iloc[pas_force_idx]

    # Pull off the data
    twitch_data = dict()
    
    max_pas_force_idx = np.argmax(pas_force_array)
    twitch_data['pas_force'] = float(pas_force_array[max_pas_force_idx])
    twitch_data['pas_force_time_s'] = float(pas_force_time_s_array[max_pas_force_idx])

    twitch_data['max_force'] = float(np.max(max_force_array))
    max_force_idx = np.argmax(max_force_array)
    twitch_data['time_s_at_max_force'] = float(max_force_time_s_array[max_force_idx])

    max_thin_act_idx = np.argmax(max_thin_act_array)
    twitch_data['time_s_at_max_thin_act'] = float(max_thin_act_time_s_array[max_thin_act_idx])

    # print(twitch_data)
    # exit(1)

    return twitch_data
