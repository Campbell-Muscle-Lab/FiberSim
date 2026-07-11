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

def twitch(json_analysis_file_string, char_index = 0):
    """ Runs a twitch protocol """

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

    # Check for protocols, which are needed for twitches
    # If not present, break out
    if not ('protocol' in char_dict):
        print('\nError: \'Protocol\' not found in characterization[%i] in %s' % 
              (char_index, json_analysis_file_string))
        print('Now exiting')
        exit(1)

    # Make the protocols
    protocol_file_strings = prepare_protocols(json_analysis_file_string,
                                              char_index = char_index)

    no_of_protocols = len(protocol_file_strings)

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
                create_sim_input_and_output_dirs(top_sim_dir,
                                                    dir_counter = dir_counter)

            # Set the new_model_file_string
            new_model_file_string = str(Path(os.path.join(sim_input_dir, Path(mfs).name)).
                                        absolute().resolve())

            # Update and write it
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
                repeat_jobs = prepare_repeats(json_analysis_file_string,
                                                sim_input_dir,
                                                sim_output_dir,
                                                new_model_file_string,
                                                protocol_file_strings[prot_index],
                                                char_index = char_index,
                                                prot_index = prot_index,
                                                impose_afterloads = False)

                # Append the jobs
                for j in repeat_jobs:
                    all_jobs.append(j)

    # Create a batch

    # We will need the folder above the sample condition
    output_dir_for_batch_figs = os.path.join(top_sim_dir, 'sim_output')

    isometric_batch = dict()
    isometric_batch['FiberSim_batch'] = dict()
    isometric_batch['FiberSim_batch']['FiberCpp_exe'] = \
        return_FiberCpp_exe_dict(json_analysis_file_string)
    isometric_batch['FiberSim_batch']['job'] = all_jobs
    isometric_batch['FiberSim_batch']['batch_figures'] = \
        return_batch_figs_dict(json_analysis_file_string,
                                output_dir_for_batch_figs,
                                char_index)

    # Write the batch
    batch_file_string = os.path.join(top_sim_dir, 'batch.json')

    with open(batch_file_string, 'w') as f:
        json.dump(isometric_batch, f, indent = 4)

    # Now work out how to run it
    figures_only = False
    if (run_mode == 'figures_only'):
        figures_only = True

    batch.run_batch(batch_file_string,
                    figures_only = figures_only)