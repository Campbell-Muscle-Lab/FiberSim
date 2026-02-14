
import os
import json

import shutil

from pathlib import Path

from .characterize_utilities import \
        create_sim_input_and_output_dirs, \
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
        update_and_write_model_file

from ..batch import batch

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

    # # Get the options_dict
    # with open(orig_options_file_string, 'r') as f:
    #     options_dict = json.load(f)

    # Prep isometric and loaded sim_dirs
    isometric_sim_dir = os.path.join(top_sim_dir, 'isometric')
    prepare_simulation_dir(isometric_sim_dir, run_mode)

    loaded_sim_dir = os.path.join(top_sim_dir, 'loaded')
    prepare_simulation_dir(loaded_sim_dir, run_mode)

    # Loop through the simulations twice, the first time
    # we are in isometric mode to establish the peak time and force
    # the second time, we are in loaded mode

    for mode_i in range(2):

        # Set some stuff up
        if (mode_i == 0):
            sim_dir = isometric_sim_dir
            impose_afterloads = False

        else:
            sim_dir = loaded_sim_dir
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
                    repeat_jobs = prepare_repeats(json_analysis_file_string,
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
        batch_file_string = os.path.join(isometric_sim_dir, 'batch.json')

        with open(batch_file_string, 'w') as f:
            json.dump(isometric_batch, f, indent = 4)

        # Now run it
        batch.run_batch(batch_file_string)

        # If we are in the first loop, we have to pull off the
        # peak times and forces to use for the second loop
        if (mode_i == 0):
            # Run the analysis


    