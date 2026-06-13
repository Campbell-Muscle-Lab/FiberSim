import os
import json
import copy

from pathlib import Path
import shutil

from ..protocols import protocols as prot

def return_run_mode(json_analysis_file_string, char_index = 0):
    """ Return the run_mode for the characterization with
    a zero-based index of char_id (default is 0) """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)
    # Pull off the char_dict
    char_dict = json_dict['FiberSim_setup']['characterization'][char_index]

    # Check
    run_mode = 'default'
    if ('run_mode' in char_dict):
        run_mode = char_dict['run_mode']
        
    # Return
    return run_mode   

def return_FiberCpp_exe_dict(json_analysis_file_string):
    """ Returns a dictionary for the FiberCpp executable information """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)
    
    # Extract the FiberCpp_exe struct
    FiberCpp_exe_dict = json_dict['FiberSim_setup']['FiberCpp_exe']

    print('ff')
    print(json_analysis_file_string)
    print(json.dumps(FiberCpp_exe_dict, indent=4))


    # If we are in a relative mode, adapt the FiberCpp_exe for absolute paths
    # because the new simulations will be run from a different folder

    if (FiberCpp_exe_dict['relative_to'] == "this_file"):
        base_dir = Path(json_analysis_file_string).parent.absolute()
        FiberCpp_exe_dict['relative_to'] = 'False'
        FiberCpp_exe_dict['exe_file'] = \
            str(Path(os.path.join(base_dir, FiberCpp_exe_dict['exe_file'])).resolve())

    return FiberCpp_exe_dict

def return_base_dir(json_analysis_file_string,
                    dict_key,
                    append_key = [],
                    append_two_keys = [],
                    append_string = [],
                    dict_index = 0):
    """ Return the base directory for the specified dict
        If the dict is an array, use the dict_id'th member """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    # Pull off the target_dict
    target_dict = json_dict['FiberSim_setup'][dict_key]

    if (isinstance(target_dict, list)):
        target_dict = target_dict[dict_index]

    # See if we need to adjust the filesystem paths
    if ('relative_to' in target_dict):
        if (target_dict['relative_to'] == "this_file"):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = target_dict['relative_to']
        if not (append_key == []):
            base_dir = os.path.join(base_dir, target_dict[append_key])
        if not (append_two_keys == []):
            base_dir = os.path.join(base_dir, target_dict[append_two_keys[0]][append_two_keys[1]])
        if not (append_string == []):
            base_dir = os.path.join(base_dir, append_string)
    else:
        base_dir = target_dict[append_key]
        if not (append_string == []):
            base_dir = os.path.join(base_dir, append_string)

    # Return
    return base_dir

def return_model_file_strings(json_analysis_file_string):
    """ Return a list of absolute paths to model files """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    model_dict = json_dict['FiberSim_setup']['model']

    # Get the base directory
    base_dir = return_base_dir(json_analysis_file_string, 'model')

    # Get the model file strings
    model_file_strings = []
    for mfs in model_dict['model_files']:
        model_file_strings.append(os.path.join(base_dir, mfs))

    # Return
    return model_file_strings

def return_options_file_string(json_analysis_file_string):
    """ Return an absolute path to the options file """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    model_dict = json_dict['FiberSim_setup']['model']

    # Get the base directory
    base_dir = return_base_dir(json_analysis_file_string, 'model')

    # Get the file_string
    options_file_string = os.path.join(base_dir, model_dict['options_file'])

    # Return
    return options_file_string


def return_hs_lengths(json_analysis_file_string, char_index = 0):
    """ Return the list of half-sarcomere lengths for the characterization with
    a zero-based index of char_index (default is 0) """
    
    print(json_analysis_file_string)

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)
    
     # Pull off the char_dict
    char_dict = json_dict['FiberSim_setup']['characterization'][char_index]
    
    # Check for half-sarcomere lengths in the char_dict
    # If not are specified, create a list from the model file
    if ('hs_lengths' in char_dict):
        hs_lengths = char_dict['hs_lengths']
    else:
        # We just want the hs_length in the first model, additional model files
        # are likely manipulations
        hs_lengths = []
        model_file_strings = return_model_file_strings(json_analysis_file_string);
        with open(model_file_strings[0], 'r') as f:
            model_dict = json.load(f)
             
        hs_lengths.append(model_dict['muscle']['initial_hs_length'])
        
    # Return
    return hs_lengths

def prepare_simulation_dir(sim_dir, run_mode):
    """ If the run_mode is 'figures_only', do nothing
        Otherewise, remove any existing files / folders in sim_dir
        and make sim_dir if it does not exist """

    if (run_mode == 'figures_only'):
        return

    # Prepare a clean dir
    prepare_clean_dir(sim_dir)

def prepare_clean_dir(clean_dir):
    # Clean the dir
    try:
        print('Trying to remove: %s' % clean_dir)
        shutil.rmtree(clean_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s')

    # Make the clean_dir if it does not exist
    if not os.path.exists(clean_dir):
        os.makedirs(clean_dir)


def create_sim_input_and_output_dirs(parent_dir, dir_counter = 1):
    """ Create sim_input and sim_output directories with a specified index """

    # Sim input
    sim_input_dir = os.path.join(parent_dir,
                                 'sim_input',
                                 ('%i' % dir_counter))

    if not os.path.isdir(sim_input_dir):
        os.makedirs(sim_input_dir)

    # Sim output
    sim_output_dir = os.path.join(parent_dir,
                                 'sim_output',
                                 ('%i' % dir_counter))

    if not os.path.isdir(sim_output_dir):
        os.makedirs(sim_output_dir)

    # Return
    return (sim_input_dir, sim_output_dir)

def update_and_write_model_file(old_model_file_string,
                                new_model_file_string,
                                char_dict,
                                hs_length):
    """ Update the model file with the new hs_length and
        any other parameters """

    with open(old_model_file_string, 'r') as f:
        new_mod = json.load(f)

    # Adjust hs_length
    new_mod['muscle']['initial_hs_length'] = float(hs_length)

    # Over-ride m_n if appropriate
    if ('m_n' in char_dict):
        new_mod['thick_structure']['m_n'] = char_dict['m_n']

    # Now write the new model file
    with open(new_model_file_string, 'w') as f:
        json.dump(new_mod, f, indent = 4)

def prepare_repeats(json_analysis_file_string,
                    sim_input_dir,
                    sim_output_dir,
                    model_file_string,
                    protocol_file_string = [],
                    char_index = 0,
                    prot_index = 0,
                    impose_afterloads = True):
    """ Prepare the repeats for the simulation """

    # Create an empty to hold the repeats
    repeat_jobs = []

    # Load the analysis file and pull off the char_dict
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    char_dict = json_dict['FiberSim_setup']['characterization'][char_index]

    # Pull off the options_dict
    orig_options_file_string = return_options_file_string(json_analysis_file_string)
    with open(orig_options_file_string, 'r') as f:
        orig_options_dict = json.load(f)

    # Handle randomized repeats if appropriate
    if ('randomized_repeats' in char_dict):
        rand_repeats = char_dict['randomized_repeats']
        orig_options_dict['options']['rand_seed'] = "random"
    else:
        rand_repeats = 1

    # Check for the number of loads
    try:
        no_of_loads = len(char_dict['protocol']['data'][prot_index]['afterload']['load'])
    except:
        no_of_loads = 1

    # print(json.dumps(char_dict, indent = 4))
    # exit(1)
    

    # Loop through loads
    for load_index in range(no_of_loads):

        # Loop through rand_repeats
        for rep_index in range(rand_repeats):

            # Copy the options
            rep_options_dict = copy.deepcopy(orig_options_dict)

            # Update the options file for status files
            if ('status_files' in rep_options_dict['options']):
                rep_options_dict['options']['status_files']['relative_to'] = \
                    'this_file'
                last_folder = \
                    Path(rep_options_dict['options']['status_files']['status_folder']).name
                status_folder = \
                    os.path.join(sim_output_dir,
                                    ('%s_%i_%i_r%i' % (last_folder,
                                                    prot_index + 1,
                                                    load_index + 1,
                                                    rep_index + 1)))
                # Make the folder if it does not exist
                if not os.path.isdir(status_folder):
                    os.makedirs(status_folder)

                # Update the dict
                rep_options_dict['options']['status_files']['status_folder'] = status_folder

            # Update the options files for rates
            if ( (prot_index == 0) and (load_index == 0) and (rep_index == 0) ):
                # This is the first trial so update the options to dump rates
                rep_options_dict['options']['rate_files'] = dict()
                rep_options_dict['options']['rate_files']['relative_to']= 'false'
                rep_options_dict['options']['rate_files']['file'] = \
                    str(Path(os.path.join(sim_output_dir, 'rates.json')).resolve())

            # If appropriate, add in afterload options
            if (impose_afterloads):
                if ('afterload' in char_dict['protocol']['data'][prot_index]):
                    after_dict = char_dict['protocol']['data'][prot_index]['afterload']
                    after_keys = after_dict.keys()
                    rep_options_dict['options']['afterload'] = dict()
                    for ak in after_keys:
                        rep_options_dict['options']['afterload'][ak] = \
                            after_dict[ak][load_index]

            # Create the options file string
            new_options_file = os.path.join(
                Path(sim_input_dir),
                ('sim_options_%i_%i_r%i.json' % (prot_index + 1,
                                                    load_index + 1,
                                                    rep_index + 1)))

            # Tidy it
            new_options_file = str(Path(new_options_file).resolve())

            # Write it
            with open(new_options_file, 'w') as f:
                json.dump(rep_options_dict, f, indent = 4)

            # Set the results file string
            results_file_string = os.path.join(sim_output_dir,
                                        ('sim_prot_%i_%i_r%i.txt' % (
                                                    prot_index + 1,
                                                    load_index + 1,
                                                    rep_index + 1)))

            # Tidy it
            results_file_string = str(Path(results_file_string).resolve())

            # Now we can create the job
            j = dict()
            j['relative_to'] = 'False'
            if (protocol_file_string):
                j['protocol_file'] = protocol_file_string
            j['results_file'] = results_file_string
            j['model_file'] = model_file_string
            j['options_file'] = new_options_file

            # Add the job to the list
            repeat_jobs.append(j)

    # Return
    return repeat_jobs

def prepare_protocols(json_analysis_file_string,
                      char_index = 0):

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    char_dict = json_dict['FiberSim_setup']['characterization'][char_index]

    # Do we already have protocol files
    if ('protocol_files' in char_dict['protocol']):
        prot_files = char_dict['protocol']['protocol_files']

    else:
        # We need to make them from the data

        # Set and create a directory for the protocols
        prot_dir = return_base_dir(json_analysis_file_string,
                                   'characterization',
                                   dict_index = char_index,
                                   append_string = char_dict['protocol']['protocol_folder'])

        prepare_clean_dir(prot_dir)

        # Now loop through the protocol data and create the protocol files, storing them
        # as we go
        prot_files = [];

        for (prot_i, prot_data) in enumerate(char_dict['protocol']['data']):

            # Prepare and tidy the protocol file string
            prot_file_string = os.path.join(prot_dir,
                                            ('protocol_%i.json' % (prot_i + 1)))

            prot_file_string = str(Path(prot_file_string).resolve())

            p = prot.create_twitch_protocol(
                    time_step = prot_data['time_step_s'],
                    n_points = prot_data['n_points'],
                    stimulus_times_s = prot_data['stimulus_times_s'],
                    Ca_content = prot_data['Ca_content'],
                    stimulus_duration_s = prot_data['stimulus_duration_s'],
                    k_leak = prot_data['k_leak'],
                    k_act = prot_data['k_act'],
                    k_serca = prot_data['k_serca'])

            prot.write_protocol_to_file(p, prot_file_string)

            prot_files.append(prot_file_string)

    # Return
    return prot_files

def return_batch_figs_dict(json_analysis_file_string,
                           sim_output_dir,
                           char_index = 0):
    """ Returns a dict for batch figs """

    # Load the analysis file and pull off the char_dict
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    char_dict = json_dict['FiberSim_setup']['characterization'][char_index]

    # Create a batch_figs dict
    batch_figs_dict = dict()

    # Append a rates fig
    batch_figs_dict['rates'] = []
    batch_figs_dict['rates'].append(
        return_batch_rates_fig_dict(sim_output_dir,
                                    char_dict))

    # And a superposed traces fig
    batch_figs_dict['superposed_traces'] = []
    batch_figs_dict['superposed_traces'].append(
        return_superposed_traces_fig_dict(sim_output_dir,
                                          char_dict))

    # Return
    return batch_figs_dict


def return_batch_rates_fig_dict(sim_output_dir, char_dict):
    """ Return a rates fig dict """

    fig = dict()
    fig['relative_to'] = 'False'
    fig['results_folder'] = sim_output_dir
    fig['output_image_file'] = os.path.join(sim_output_dir, 'rates')
    fig['output_image_formats'] = char_dict['output_image_formats']
    if ('formatting' in char_dict):
            fig['formatting'] = char_dict['formatting']

    return fig

def return_superposed_traces_fig_dict(sim_output_dir, char_dict):
    """ Return a superposed traces fig dict """

    fig = dict()
    fig['relative_to'] = 'False'
    fig['results_folder'] = sim_output_dir
    fig['output_image_file'] = str(Path(os.path.join(sim_output_dir, 'superposed_traces')).resolve())
    fig['output_image_formats'] = char_dict['output_image_formats']
    if ('superposed_x_ticks' in char_dict):
           fig['superposed_x_ticks'] = char_dict['superposed_x_ticks']
    if ('formatting' in char_dict):
           fig['formatting'] = char_dict['formatting']

    return fig

def update_setup_file_string_with_new_char(
    json_dict, char_dict, char_index, new_setup_file_string):
    """ Inserts the char_dict into the setup and writes the file """

    json_dict['FiberSim_setup']['characterization'][char_index] = char_dict

    with open(new_setup_file_string, 'w') as f:
        json.dump(json_dict, f, indent=4)
