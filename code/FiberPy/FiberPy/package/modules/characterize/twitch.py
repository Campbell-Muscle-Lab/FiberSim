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
                    create_sim_input_and_output_dirs(sim_dir,
                                                     dir_counter = dir_counter)
