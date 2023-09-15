# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 12:47:26 2022

@author: kscamp3
"""
import os
import json
import shutil
import copy
import re


from natsort import natsorted

import numpy as np
import pandas as pd

from pathlib import Path

from ..protocols import protocols as prot
from ..batch import batch

# from .characterize_functions import characterize_fv_with_pCa_and_isometric_force


def characterize_model(json_analysis_file_string):
    """ Code takes a json struct that includes a model file, and run the
        analyses described in the file """
        
    # Check for the analysis file
    if (not json_analysis_file_string):
        print('deduce_fv_properties: no analysis file specified')
        exit(1)
    
    # Load it
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        anal_struct = json_data['FiberSim_characterization']
        
    # If there is a manipulations section, use that to create the
    # appropriate models
    if ("manipulations" in anal_struct['model']):
        json_analysis_file_string = \
            generate_model_files(json_analysis_file_string)
            
    print(json_analysis_file_string)
   
    # Pull off the characterization tasks
    char_struct = anal_struct['characterization']
    for ch in char_struct:
        if (ch['type'] == 'pCa_length_control'):
            deduce_pCa_length_control_properties(json_analysis_file_string,
                                                 pCa_struct = ch)
        if (ch['type'] == 'force_velocity'):
            deduce_fv_properties(json_analysis_file_string,
                                 fv_struct =ch)
        
        if (ch['type'] == 'freeform'):
            deduce_freeform_properties(json_analysis_file_string,
                                       freeform_struct = ch)
        
        if (ch['type'] == 'fv_with_pCa_and_isometric_force'):
            characterize_fv_with_pCa_and_isometric_force(
                json_analysis_file_string,
                ch)
        
def generate_model_files(json_analysis_file_string):
    """ Clones base model with modifications to facilitate comparisons """
    
    # First load the file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        model_struct = json_data['FiberSim_characterization']['model']
    
    # Deduce the base model file string
    if (model_struct['relative_to'] == 'this_file'):
        model_working_dir = Path(json_analysis_file_string).parent.absolute()
        base_model_file_string = os.path.join(model_working_dir,
                                              model_struct['manipulations']['base_model'])
        base_model_file_string = str(Path(base_model_file_string).resolve())
    else:
        base_model_file_string = model_struct['manipulations']['base_model']

    # Now deduce where to put the adjusted model files
    if (model_struct['relative_to'] == 'this_file'):
        generated_dir = os.path.join(Path(json_analysis_file_string).parent.absolute(),
                                     model_struct['manipulations']['generated_folder'])
        generated_dir = str(Path(generated_dir).resolve())
    else:
        generated_dir = model_struct['manipulations']['generated_folder']
        
    # Clean the generated dir
    try:
        print('Trying to remove %s' % generated_dir)
        shutil.rmtree(generated_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s' % (generated_dir, e.strerror))
        
    if not os.path.isdir(generated_dir):
        os.makedirs(generated_dir)
        
    # Now copy the sim_options file across
    if (model_struct['relative_to'] == 'this_file'):
        temp_dir = Path(json_analysis_file_string).parent.absolute()
        orig_options_file = os.path.join(temp_dir,
                                         model_struct['options_file'])
        new_options_file = os.path.join(generated_dir, model_struct['options_file'])
    else:
        orig_options_file = model_struct['options_file']
        print('Needs more work on sim_options copying')
        exit(1)
    
    shutil.copy(orig_options_file, new_options_file)

    # Load up the base model
    with open(base_model_file_string, 'r') as f:
        base_model = json.load(f)

    # Now work out how many adjustments you need to make
    adjustments = model_struct['manipulations']['adjustments']
    no_of_models = len(adjustments[0]['multipliers'])
    
    generated_models = []
        
    # Loop through them
    for i in range(no_of_models):
        
        # Copy the base model
        adj_model = copy.deepcopy(base_model)
                
        for (j,a) in enumerate(adjustments):
            
            value = a['base_value'] * a['multipliers'][i]
            
            if (a['output_type'] == 'int'):
                adj_model[a['class']][a['variable']] = int(value)
                
            if (a['output_type'] == 'float'):
                adj_model[a['class']][a['variable']] = float(value)
                
            # Check for NaN
            if (np.isnan(value)):
                adj_model[a['class']][a['variable']] = 'null'
        
    
        # Now generate the model file string
        model_file_string = 'model_%i.json' % (i+1)
        
        # We need the full path to write it to disk
        adj_model_file_string = os.path.join(generated_dir,
                                             model_file_string)

        with open(adj_model_file_string, 'w') as f:
            json.dump(adj_model, f, indent=4)

        # Append the model files
        generated_models.append(model_file_string)
        
    # Update the set up file
    
    # Add in the model files
    json_data['FiberSim_characterization']['model']['model_files'] = generated_models
    
    # Delete the adjustments
    del(json_data['FiberSim_characterization']['model']['manipulations'])
    
    # Generate a new setup file string
    generated_setup_file_string = os.path.join(generated_dir,
                                               'generated_setup.json')
    
    # Write to file
    with open(generated_setup_file_string, 'w') as f:
        json.dump(json_data, f, indent=4)
        
    # Return the new filename
    return (generated_setup_file_string)
            
            
def deduce_pCa_length_control_properties(json_analysis_file_string,
                                         pCa_struct = []):
    """ Code runs pCa analysis """

    # Potentially switch off simulations
    figures_only = False
    if ('figures_only' in pCa_struct):
        if (pCa_struct['figures_only'] == "True"):
            figures_only = True

    trace_figures_on = True
    if ('trace_figures_on' in pCa_struct):
        if (pCa_struct['trace_figures_on'] == 'False'):
            trace_figures_on = False

    # Load the file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        char_struct = json_data['FiberSim_characterization']
    
    # Pull off the components
    FiberCpp_exe_struct = char_struct['FiberCpp_exe']
    model_struct = char_struct['model']
    
    # Create an isometric batch to run the isometric test
    pCa_lc_b = dict()
    
    # Turn the FiberCpp_exe into absolute paths because the new instruction
    # file will be in a different place, add to new batch
    cpp_exe = dict()
    if ('relative_to' in FiberCpp_exe_struct):
        if (FiberCpp_exe_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = FiberCpp_exe_struct['relative_to']
        cpp_exe['exe_file'] = os.path.join(base_dir,
                                           FiberCpp_exe_struct['exe_file'])
    else:
        cpp_exe['exe_file'] = FiberCpp_exe_struct['exe_file']

    pCa_lc_b['FiberCpp_exe'] = cpp_exe

    pCa_lc_b['job'] = []
    
    # Check for half-sarcomere lengths in the pCa_struct
    # If none are specified, create an hsl array from the model file
    if ('hs_lengths' in pCa_struct):
        hs_lengths = pCa_struct['hs_lengths']
    else:
        model_file = os.path.join(Path(json_analysis_file_string).parent.absolute(),
                                  model_struct['model_files'][0])
        with open(model_file, 'r') as f:
            m = json.load(f)
            hs_lengths = np.array([m['muscle']['initial_hs_length']])
            
    if ('length_step_nm' in pCa_struct):
        length_steps = pCa_struct['length_step_nm']
    else:
        length_steps = [0]
    
    # Deduce the base_dir
    if ('relative_to' in pCa_struct):
        if (pCa_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = pCa_struct['relative_to']
        base_dir = os.path.join(base_dir, pCa_struct['sim_folder'])
    else:
        base_dir = pCa_struct['sim_folder']        
    
    # If you are running simulations, delete the existing structure
    if not figures_only:
        try:
            print('Trying to remove %s' % base_dir)
            shutil.rmtree(base_dir, ignore_errors=True)
        except OSError as e:
            print("Error: %s : %s" % (base_dir, e.strerror))
    
    # Set up dir_counter
    dir_counter = 0
    
    # Loop through the model files
    for i, mod_f in enumerate(model_struct['model_files']):
        
        # Now loop through the half-sarcomere lengths
        for j, hsl in enumerate(hs_lengths):
            
            # Now loop through the length_steps
            for k, length_step in enumerate(length_steps):
                
                # Update dir_counter
                dir_counter = dir_counter + 1
                
                # Create the input and output directories           
                sim_input_dir = os.path.join(base_dir,
                                             'sim_input',
                                             ('%i' % dir_counter))
                if not os.path.isdir(sim_input_dir):
                    os.makedirs(sim_input_dir)
                    
                sim_output_dir = os.path.join(base_dir,
                                              'sim_output',
                                              ('%i' % dir_counter))
                if not os.path.isdir(sim_output_dir):
                    os.makedirs(sim_output_dir)
                   
                # Copy the model and options files to the sim_input dir
                # adjusting half-sarcomere length as appropriate
                if (model_struct['relative_to'] == 'this_file'):
                    model_dir = Path(json_analysis_file_string).parent.absolute()
                    orig_model_file = os.path.join(model_dir, mod_f)
                else:
                    orig_model_file = mod_f
                
                # Adjust hsl by loading model, adjusting hsl and re-writing
                with open(orig_model_file, 'r') as f:
                    m = json.load(f)
                    m['muscle']['initial_hs_length'] = float(hsl)
                    
                    # Over-ride m_n if appropriate
                    if ('m_n' in pCa_struct):
                        m['thick_structure']['m_n'] = pCa_struct['m_n']
                        
                fn = re.split('/|\\\\', orig_model_file)[-1]
                model_file = os.path.join(sim_input_dir, fn)
    
                with open(model_file, 'w') as f:
                    json.dump(m, f, indent=4)
                
                # Work out the path for the base options file
                if (model_struct['relative_to'] == 'this_file'):
                    model_dir = Path(json_analysis_file_string).parent.absolute()
                    orig_options_file = os.path.join(model_dir,
                                                     model_struct['options_file'])
                else:
                    orig_options_file = model_struct['options_file']
    
                # Load the options data
                with open(orig_options_file, 'r') as f:
                    orig_options_data = json.load(f)
    
                # Adjust the options struct if we have randomized repeats
                if ('randomized_repeats' in pCa_struct):
                    rand_repeats = pCa_struct['randomized_repeats']
                    orig_options_data['options']['rand_seed'] = "random"
                else:
                    rand_repeats = 1
               
                # Loop through the pCa values
                for pCa_counter,pCa in enumerate(pCa_struct['pCa_values']):
                    
                    # Loop through the rand_repeats, creating a job for each repeat
                    for rep in range(rand_repeats):
                        
                        # Copy the orig_options_struct for local changes
                        # within the rep
                        rep_options_data = copy.deepcopy(orig_options_data)
                        
                        # Update the options file to dump to a local directory
                        if ('status_files' in rep_options_data['options']):
                            rep_options_data['options']['status_files']['status_folder'] = \
                                os.path.join(sim_output_dir,
                                            ('%s_%i_r%i' % (rep_options_data['options']['status_files']['status_folder'],
                                                          (pCa_counter + 1), (rep+1))))
                                
                        if ((pCa_counter==0) and (rep == 0)):
                            # If it is the first pCa_value for the model and the first rep,
                            # update the options file to dump rates
                            rep_options_data['options']['rate_files'] = dict()
                            rep_options_data['options']['rate_files']['relative_to'] = \
                                'false'
                            rep_options_data['options']['rate_files']['file'] = \
                                os.path.join(sim_output_dir, 'rates.json')
                        
                        # Create the new options file
                        options_file = os.path.join(
                            sim_input_dir,
                            ('sim_options_%i_r%i.json' % (pCa_counter+1,
                                                          rep+1)))
            
                        with open(options_file, 'w') as f:
                                json.dump(rep_options_data, f, indent=4)
        
                        # Make a protocol, thinking about whether we need smaller
                        # time steps for k_tr
                        if not ('k_tr_start_s' in pCa_struct):
                            n_points = int(pCa_struct['sim_duration_s'] /
                                           pCa_struct['time_step_s'])
                            dt = pCa_struct['time_step_s'] * np.ones(n_points)
                            pCa_vector = pCa * np.ones(n_points)
                            delta_hsl = np.zeros(n_points)
                            mode_vector = -2 * np.ones(n_points)
                        else:
                            # Pre-phase
                            if not ('length_step_nm' in pCa_struct):
                                pre_points = int(pCa_struct['k_tr_start_s'] /
                                                 pCa_struct['time_step_s'])
                                pre_dt = pCa_struct['time_step_s'] * np.ones(pre_points)
                                pre_pCa = pCa * np.ones(pre_points)
                                pre_delta_hsl = np.zeros(pre_points)
                                pre_mode_vector = -2 * np.ones(pre_points)
                            else:
                                # Break up pre_ktr period into before length step,
                                # during length step and after length step
                                pre_1_points = int(pCa_struct['length_step_s'] /
                                                   pCa_struct['time_step_s'])
                                pre_1_dt = pCa_struct['time_step_s'] * np.ones(pre_1_points)
                                pre_1_delta_hsl = np.zeros(pre_1_points)
                                pre_1_mode_vector = -2 * np.ones(pre_1_points)
                                
                                ramp_time_step = pCa_struct['time_step_s']  / 10
                                step_ramp_points = int(pCa_struct['length_step_ramp_s'] /
                                                       ramp_time_step)
                                step_ramp_dt = ramp_time_step * np.ones(step_ramp_points)
                                step_ramp_inc = length_step / float(step_ramp_points)
                                step_ramp_delta_hsl = step_ramp_inc * np.ones(step_ramp_points)
                                step_ramp_mode_vector = -1 * np.ones(step_ramp_points)
                                
                                pre_2_points = int((pCa_struct['k_tr_start_s'] -
                                                        pCa_struct['length_step_s'] -
                                                        pCa_struct['length_step_ramp_s']) /
                                                   pCa_struct['time_step_s'])
                                pre_2_dt = pCa_struct['time_step_s'] * np.ones(pre_2_points)
                                pre_2_delta_hsl = np.zeros(pre_2_points)
                                pre_2_mode_vector = -2 * np.ones(pre_2_points)
                                
                                # Assemble into pre arrays
                                pre_points = pre_1_points + step_ramp_points + pre_2_points
                                pre_dt = np.hstack((pre_1_dt,
                                                    step_ramp_dt,
                                                    pre_2_dt))
                                pre_pCa = pCa * np.ones(pre_points)
                                pre_delta_hsl = np.hstack((pre_1_delta_hsl,
                                                           step_ramp_delta_hsl,
                                                           pre_2_delta_hsl))
                                pre_mode_vector = np.hstack((pre_1_mode_vector,
                                                             step_ramp_mode_vector,
                                                             pre_2_mode_vector))
                                
                            
                            # k_tr
                            k_tr_time_step = pCa_struct['time_step_s'] / 10
                            k_tr_points = int(pCa_struct['k_tr_duration_s'] / k_tr_time_step)
                            k_tr_ramp_points = int(pCa_struct['k_tr_ramp_s'] / k_tr_time_step)
                            ramp_inc = pCa_struct['k_tr_magnitude_nm'] / \
                                        float(k_tr_ramp_points)
                            
                            k_tr_dt = k_tr_time_step * np.ones(k_tr_points)
                            k_tr_pCa = pCa * np.ones(k_tr_points)
                            k_tr_delta_hsl = np.zeros(k_tr_points)
                            vi = np.arange(0, k_tr_ramp_points+1)
                            k_tr_delta_hsl[vi] = -ramp_inc
                            vi = np.arange(k_tr_points-1-k_tr_ramp_points, k_tr_points)
                            k_tr_delta_hsl[vi] = ramp_inc
                            k_tr_mode_vector = -np.ones(k_tr_points)
    
                            # Post
                            post_points = int((pCa_struct['sim_duration_s'] -
                                               pCa_struct['k_tr_start_s'] -
                                               pCa_struct['k_tr_duration_s']) /
                                              pCa_struct['time_step_s'])
                            post_dt = pCa_struct['time_step_s'] * np.ones(post_points)
                            post_pCa = pCa * np.ones(post_points)
                            post_delta_hsl = np.zeros(post_points)
                            post_mode_vector = -2 * np.ones(post_points)
                            
                            # Stack together
                            dt = np.hstack((pre_dt, k_tr_dt, post_dt))
                            pCa_vector = np.hstack((pre_pCa, k_tr_pCa, post_pCa))
                            delta_hsl = np.hstack((pre_delta_hsl,
                                                   k_tr_delta_hsl,
                                                   post_delta_hsl))
                            mode_vector = np.hstack((pre_mode_vector,
                                                     k_tr_mode_vector,
                                                     post_mode_vector))
                            
                        # Now adjust the pCa_vector for step_up and step_down
                        # if required
                        t = np.cumsum(dt)
                        if ('pCa_start' not in pCa_struct):
                            pCa_start = 9.0
                        else:
                            pCa_start = pCa_struct['pCa_start']
                        
                        if ('pCa_stop' not in pCa_struct):
                            pCa_stop = 9.0
                        else:
                            pCa_stop = pCa_struct['pCa_stop']
                            
                        if ('pCa_step_up_s' in pCa_struct):
                            pCa_vector[t < pCa_struct['pCa_step_up_s']] = pCa_start
                        
                        if ('pCa_step_down_s' in pCa_struct):
                            pCa_vector[t > pCa_struct['pCa_step_down_s']] = pCa_stop
                        
                        # Now make the protocol
                        df = pd.DataFrame({'dt': dt,
                                           'pCa': pCa_vector,
                                           'delta_hsl': delta_hsl,
                                           'mode': mode_vector})
                        
                        prot_file_string = os.path.join(sim_input_dir,
                                                        ('prot_pCa_%.0f_s_%.0f_r%i.txt' %
                                                         (10*pCa, 10*length_step, rep+1)))
                        
                        # Write protocol if required
                        if not figures_only:
                            prot.write_protocol_to_file(df, prot_file_string)
                    
                        # Create the job
                        j = dict()
                        j['relative_to'] = 'False'
                        j['protocol_file'] = prot_file_string
                        j['results_file'] = os.path.join(sim_output_dir,
                                                         ('sim_pCa_%.0f_s_%.0f_r%i.txt' %
                                                          (10*pCa, 10*length_step, rep+1)))
                        j['model_file'] = model_file
                        j['options_file'] = options_file
                    
                        # If required, create an output_handler and add it to
                        # the job
                        if (trace_figures_on == True):
                            print('Need to fix trace figures')
                            exit(1)
                            # Create the structure for the output handler
                            oh = dict()
                            oh['templated_images'] = []
                            tf = dict()
                            tf['relative_to'] = 'this_file'
                            
                            tf['template_file_string'] = os.path.join(
                                                            '..',
                                                            base_dir,
                                                            'template',
                                                            'template_summary.json')
                            tf['output_file_string'] = os.path.join(
                                                            sim_output_dir,
                                                            ('sim_pCa_%.0f_r%i' %
                                                             (10*pCa, rep+1)))
                            tf['output_image_formats'] = pCa_struct['output_image_formats']
                            oh['templated_images'].append(tf)
                            
                            # Now add it to the job, and write it to file
                            j['output_handler_file'] = os.path.join(
                                                        sim_input_dir,
                                                        ('output_handler_pCa_%.0f_r%i.json' %
                                                            (10*pCa, rep+1)))
                            
                            with open(j['output_handler_file'], 'w') as f:
                                json.dump(oh, f, indent=4)        
            
                        pCa_lc_b['job'].append(j)

    # Now create the analysis section
    batch_figs = dict()
    
    # Deduce the output dir
    output_dir = str(Path(sim_output_dir).parent)
    
    # Create a dict for functino return values
    func_output = dict()

    # pCa curves
    batch_figs['pCa_curves'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = output_dir
    fig['data_field'] = 'hs_1_force'
    fig['output_data_file_string'] = os.path.join(output_dir,
                                                  'pCa_analysis.xlsx')
    fig['output_image_file'] = os.path.join(output_dir,
                                            'force_pCa')
    fig['output_image_formats'] = pCa_struct['output_image_formats']

    if ('formatting' in pCa_struct):
        fig['formatting'] = pCa_struct['formatting']
    else:
        fig['formatting'] = dict()
        fig['formatting']['y_axis_label'] = 'Force (N m$^{\\mathregular{-2}}$)'
    batch_figs['pCa_curves'].append(fig)

    func_output['pCa_analysis_file_string'] = fig['output_data_file_string']

    # # pCa curves - NORMALIZED
    # fig = dict()
    # fig['relative_to'] = "False"
    # fig['results_folder'] = os.path.join(output_dir)
    # fig['data_field'] = 'force'
    # fig['output_image_file'] = os.path.join(output_dir,
    #                                         'force_pCa_normalized')
    # fig['output_image_formats'] = pCa_struct['output_image_formats']

    # fig['formatting'] = dict()
    # fig['formatting']['y_axis_label'] = 'Normalized \n force'
    # fig['formatting']['y_normalized'] = 'True'
    # fig['formatting']['y_label_pad'] = 20

    # batch_figs['pCa_curves'].append(fig)

    # Rates
    batch_figs['rates'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = output_dir
    fig['output_image_file'] = os.path.join(output_dir, 'rates')
    fig['output_image_formats'] = pCa_struct['output_image_formats']
    
    if ('formatting' in pCa_struct):
        fig['formatting'] = pCa_struct['formatting']
    
    batch_figs['rates'].append(fig)

    # Superposed traces
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = output_dir
    fig['output_image_file'] = os.path.join(output_dir,
                                            'superposed_traces')
    fig['output_image_formats'] = pCa_struct['output_image_formats']
    if ('superposed_x_ticks' in pCa_struct):
        fig['superposed_x_ticks'] = pCa_struct['superposed_x_ticks']
    if ('formatting' in pCa_struct):
        fig['formatting'] = pCa_struct['formatting']        
    batch_figs['superposed_traces'].append(fig)
    
    pCa_lc_b['batch_figures'] = batch_figs

    # k_tr
    if ('k_tr_start_s' in pCa_struct):
        batch_figs['k_tr_analysis'] = []
        fig = dict()
        fig['relative_to'] = "False"
        fig['results_folder'] = output_dir
        fig['output_data_file_string'] = os.path.join(output_dir,
                                                      'k_tr_analysis.xlsx')
        fig['output_image_file'] = os.path.join(output_dir,
                                                'k_tr_analysis')
        fig['k_tr_fit_time_s'] = pCa_struct['k_tr_fit_time_s']
        fig['output_image_formats'] = pCa_struct['output_image_formats']
        if ('k_tr_ticks' in pCa_struct):
            fig['k_tr_ticks'] = pCa_struct['k_tr_ticks']
        if ('formatting' in pCa_struct):
            fig['formatting'] = pCa_struct['formatting']
        batch_figs['k_tr_analysis'].append(fig)

    # Now insert isometric_b into a full batch structure
    pCa_lc_batch = dict()
    pCa_lc_batch['FiberSim_batch'] = pCa_lc_b

    # Write the batch to file
    base_dir = str(Path(output_dir).parent)
    pCa_lc_batch_file = os.path.join(base_dir,
                                     'batch_pCa.json')
    print(pCa_lc_batch_file)
    with open(pCa_lc_batch_file, 'w') as f:
        json.dump(pCa_lc_batch, f, indent=4)
        
    # Now run the isometric batch
    batch.run_batch(pCa_lc_batch_file, figures_only=figures_only)
    
    # Add the batch to the output
    func_output['sim_batch'] = pCa_lc_batch    
    
    # Return the output
    return func_output
    

def deduce_fv_properties(json_analysis_file_string,
                          fv_struct = []):
    """ Code runs force-velocity analysis """

    # Potentially switch off simulations
    figures_only = False
    if ('figures_only' in fv_struct):
        if (fv_struct['figures_only'] == "True"):
            figures_only = True

    trace_figures_on = True
    if ('trace_figures_on' in fv_struct):
        if (fv_struct['trace_figures_on'] == 'False'):
            trace_figures_on = False

    # Load the file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        anal_struct = json_data['FiberSim_characterization']
    
    # Pull off the components
    FiberCpp_exe_struct = anal_struct['FiberCpp_exe']
    model_struct = anal_struct['model']

    # First create a batch to run an isometric test. We use this to
    # get the force values for the fv tests

    # Create the isometric batch
    isometric_b = dict()
    
    # Turn the FiberCpp_exe into absolute paths because the new instruction file
    # will be in a different place
    if (FiberCpp_exe_struct['relative_to'] == 'this_file'):
        base_dir = Path(json_analysis_file_string).parent.absolute()
        FiberCpp_exe_struct['relative_to'] = 'False'
        FiberCpp_exe_struct['exe_file'] = \
            os.path.join(base_dir, FiberCpp_exe_struct['exe_file'])
    isometric_b['FiberCpp_exe'] = FiberCpp_exe_struct

    # Check for half-sarcomere lengths in the fv_struct
    # If none are specified, create an hsl array from the model file
    if ('hs_lengths' in fv_struct):
        hs_lengths = fv_struct['hs_lengths']
    else:
        model_file = os.path.join(Path(json_analysis_file_string).parent.absolute(),
                                  model_struct['model_files'][0])
        with open(model_file, 'r') as f:
            m = json.load(f)
            hs_lengths = np.array([m['muscle']['initial_hs_length']])
    
    
    # Deduce the base_dir
    if ('relative_to' in fv_struct):
        if (fv_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = fv_struct['relative_to']
        base_dir = os.path.join(base_dir, fv_struct['sim_folder'])
    else:
        base_dir = fv_struct['sim_folder']   
    
    # If you are running simulations, delete the existing sim folders
    if not figures_only:
        try:
            print('Trying to remove %s' % base_dir)
            shutil.rmtree(base_dir, ignore_errors=True)
        except OSError as e:
            print("Error: %s : %s" % (base_dir, e.strerror))
           
    if (fv_struct['pCa'] == 'pCa_50'):
        # Run a force-pCa curve and deduce pCa50 by repurposing the
        # deduce_pCa_length_control_properties function
        
        # Adjust the output directory
        pCa_struct = copy.deepcopy(fv_struct)
        pCa_struct['sim_folder'] = os.path.join(base_dir,
                                                'pCa_curves')
         
        pCa_output = deduce_pCa_length_control_properties(
                        json_analysis_file_string,
                        pCa_struct = pCa_struct)
        
        # Now pull off the pCa_values
        # First get the sheet names
        fi = pd.ExcelFile(pCa_output['pCa_analysis_file_string'])
        excel_sheets = fi.sheet_names
        # Restrict to curves and sort
        excel_sheets.remove('simulation_data')
        excel_sheets = natsorted(excel_sheets)
        
        force_pCa_fit = []
        
        for i,c in enumerate(excel_sheets):
            d = pd.read_excel(pCa_output['pCa_analysis_file_string'],
                              sheet_name = c)
            
            c_fit = dict()
            c_fit['pCa_50'] = d['pCa_50'][0]
            c_fit['n_H'] = d['n_H'][0]                
            c_fit['y_min'] = d['y_min'][0]
            c_fit['y_amp'] = d['y_amp'][0]
            
            # Calculate force at pCa50
            c_fit['f_at_pCa_50'] = c_fit['y_min'] + (0.5 * c_fit['y_amp'])
            
            force_pCa_fit.append(c_fit)
            
        # Set the isometric jobs
        isometric_jobs = pCa_output['sim_batch']['FiberSim_batch']['job']
        
    else:
        # Run an isometric simulation for the specified pCa value
        
        # Create an entry for the jobs
        isometric_b['job'] = []
            
        # Set up dir_counter
        dir_counter = 0
        
        # Loop through the model files
        for i, mod_f in enumerate(model_struct['model_files']):
            
            # Now loop through the half-sarcomere lengths
            for j, hsl in enumerate(hs_lengths):
                
                # Update dir counter
                dir_counter = dir_counter + 1
                
                if (model_struct['relative_to'] == 'this_file'):
                    model_dir = Path(json_analysis_file_string).parent.absolute()
                else:
                    model_dir = ''
                    
                # Create the input and output directories
                sim_input_dir = os.path.join(base_dir,
                                             'isometric',
                                             'sim_input',
                                             ('%i' % dir_counter))
                if not os.path.isdir(sim_input_dir):
                    os.makedirs(sim_input_dir)
                    
                sim_output_dir = os.path.join(base_dir,
                                             'isometric',
                                             'sim_output',
                                             ('%i' % dir_counter))
                
                # Delete any existing files in output if running simulations
                if not figures_only:
                    try:
                        print('Trying to remove %s' % sim_output_dir)
                        shutil.rmtree(sim_output_dir, ignore_errors=True)
                    except OSError as e:
                        print("Error: %s : %s" % (sim_output_dir, e.strerror))
                
                if not os.path.isdir(sim_output_dir):
                    os.makedirs(sim_output_dir)
                
                orig_model_file = os.path.join(model_dir, mod_f)
                    
                # Adjust hsl by loading model, adjusting hsl and re-writing
                with open(orig_model_file, 'r') as f:
                    m = json.load(f)
                    m['muscle']['initial_hs_length'] = float(hsl)
                    
                    # Over-ride m_n if appropriate
                    if ('m_n' in fv_struct):
                        m['thick_structure']['m_n'] = fv_struct['m_n']
                
                # fn = orig_model_file.split('/')[-1]
                fn = re.split('/|\\\\', orig_model_file)[-1]
                iso_model_file = os.path.join(sim_input_dir, fn)
                
                with open(iso_model_file, 'w') as f:
                    json.dump(m, f, indent=4)
                        
                # Copy the base options file
                orig_options_file = os.path.join(model_dir,
                                                 model_struct['options_file'])

                # Load the options
                with open(orig_options_file, 'r') as f:
                    orig_options_data = json.load(f)

                
                # Copy the orig_options for local changes within the rep
                rep_options_data = copy.deepcopy(orig_options_data)
                    
                # Update the options file to dump to a local directory
                if ('status_files' in rep_options_data['options']):
                    rep_options_data['options']['status_files']['relative_to'] = \
                        'false'
                    rep_options_data['options']['status_files']['status_folder'] = \
                        os.path.join(sim_output_dir,
                                      ('%s' % (rep_options_data['options']['status_files']['status_folder'])))
                    
                    # Make the status folder if required
                    test_dir = rep_options_data['options']['status_files']['status_folder']
                    if not os.path.isdir(test_dir):
                        os.makedirs(test_dir)
               
                fn = re.split('/|\\\\', orig_options_file)[-1]
                iso_options_file = os.path.join(sim_input_dir, fn)
                
                with open(iso_options_file, 'w') as f:
                    json.dump(rep_options_data, f, indent=4)
            
                # Create a length control protocol and write it to file
                df = prot.create_length_control_protocol(
                                        time_step = fv_struct['time_step_s'],
                                        step_pCa = fv_struct['pCa'],
                                        n_points = int(fv_struct['sim_duration_s'] /
                                                            fv_struct['time_step_s']))
                prot_file_string = os.path.join(sim_input_dir,
                                                'prot_iso_pCa_%.0f.txt' % (10 * fv_struct['pCa']))
                
                # Write the protocol, unless you are running figures only
                if not figures_only:
                    prot.write_protocol_to_file(df, prot_file_string);
         
                # Now create a job
                j = dict()
                j['relative_to'] = 'False'
                j['protocol_file'] = prot_file_string
                j['results_file'] = os.path.join(sim_output_dir,
                                                 'sim_pCa_%.0f.txt' % (10 * fv_struct['pCa']))
                j['model_file'] = iso_model_file
                j['options_file'] = iso_options_file
        
                # If required, create an output_handler and add it to
                # the job
                if (trace_figures_on == True):
                    print('Need to fix paths for trace figures')
                    exit(1)
                    # Create the structure for the output handler
                    oh = dict()
                    oh['templated_images'] = []
                    tf = dict()
                    tf['relative_to'] = 'this_file'
                    tf['template_file_string'] = os.path.join(
                                                    '..',
                                                    base_dir,
                                                    'template',
                                                    'template_summary.json')
                    tf['output_file_string'] = os.path.join(
                                                    base_dir, fv_struct['sim_folder'],
                                                    'isometric', 'sim_output',
                                                    ('%i' % dir_counter),
                                                    'sim_pCa_%.0f' % (10 * fv_struct['pCa']))
                    tf['output_image_formats'] = fv_struct['output_image_formats']
                    oh['templated_images'].append(tf)
                        
                    # Now add it to the job, and write it to file
                    j['output_handler_file'] = os.path.join(
                                                sim_input_dir,
                                                'output_handler_iso_pCa_%.0f.json' %
                                                    (10 * fv_struct['pCa']))
                    
                    with open(j['output_handler_file'], 'w') as f:
                        json.dump(oh, f, indent=4)
                        
                
    
                isometric_b['job'].append(j)
    
        # Now insert iso_b into a full batch structure
        isometric_batch = dict()
        isometric_batch['FiberSim_batch'] = isometric_b    
        
        # Set the batch output dir
        # Create the folders for the analysis
        isometric_batch_file = os.path.join(base_dir,
                                            'isometric',
                                            'batch_isometric.json')

        with open(isometric_batch_file, 'w') as f:
            json.dump(isometric_batch, f, indent=4)
            
        # Now run the isometric batch
        batch.run_batch(isometric_batch_file, figures_only=figures_only)
        
        # Save the isometric jobs
        isometric_jobs = isometric_batch['FiberSim_batch']['job']

    # Switch to isotonic stuff
    # First create the isotonic batch dict
    isotonic_b = dict()
    isotonic_b['FiberCpp_exe'] = FiberCpp_exe_struct
    isotonic_b['job'] = []
    
    # Now cycle thought the isometric jobs, generating an isotonic suite
    # for each one
    
    # Set up the dir counter
    dir_counter = 0
    
    # Loop through the model files
    for i, mod_f in enumerate(model_struct['model_files']):
        
        # Now loop through the half-sarcomere lengths
        for j, hsl in enumerate(hs_lengths):

            # Update dir_counter
            dir_counter = dir_counter + 1
            
            if (fv_struct['pCa'] == "pCa_50"):
                # Deduce data from the curve fit
                c_fit = force_pCa_fit[dir_counter - 1]
                isometric_force = c_fit['f_at_pCa_50']
                isometric_job_index = (dir_counter - 1) * \
                                        len(fv_struct['pCa_values']) * \
                                        fv_struct['randomized_repeats']

            else:
                # Pull off the isoric force for the preceding job
                isometric_job_index = dir_counter - 1
                results_file_string = isometric_jobs[isometric_job_index]['results_file']
                sim_data = pd.read_csv(results_file_string, sep='\t')
                 # take the mean force over last 50 points
                isometric_force = sim_data['force'].iloc[-50:-1].mean()
    
            # Create folders for the isotonic sim_input and sim_output
            sim_input_dir = os.path.join(base_dir,
                                         'isotonic',
                                         'sim_input',
                                         ('%i' % dir_counter))
            if not os.path.isdir(sim_input_dir):
                os.makedirs(sim_input_dir)
                
            sim_output_dir = os.path.join(base_dir,
                                          'isotonic',
                                          'sim_output',
                                          ('%i' % dir_counter))
            
            # Delete any existing files in output if running simulations
            if not figures_only:
                try:
                    print('Trying to remove %s' % sim_output_dir)
                    shutil.rmtree(sim_output_dir, ignore_errors=True)
                except OSError as e:
                    print("Error: %s : %s" % (sim_output_dir, e.strerror))
            
            if not os.path.isdir(sim_output_dir):
                os.makedirs(sim_output_dir)
       
            # Copy the model file to the sim_input dir
            orig_model_file = isometric_jobs[isometric_job_index]['model_file']
            # fn = orig_model_file.split('\\')[-1]
            fn = re.split('/|\\\\', orig_model_file)[-1]
            isotonic_model_file = os.path.join(sim_input_dir, fn)
            shutil.copyfile(orig_model_file, isotonic_model_file)
            
            # Load the options
            orig_options_file = isometric_jobs[isometric_job_index]['options_file']
            with open(orig_options_file, 'r') as f:
                orig_options_data = json.load(f)
            
            # Adjust the options struct if we have randomized repeats
            if ('randomized_repeats' in fv_struct):
                rand_repeats = fv_struct['randomized_repeats']
                orig_options_data['options']['rand_seed'] = "random"
            else:
                rand_repeats = 1
    
            # Loop through the isotonic forces
            for (k, rel_f) in enumerate(fv_struct['rel_isotonic_forces']):
    
                # Loop through the rand repeats
                for rep in range(rand_repeats):
                    
                    # Copy the orig_options for local changes within the rep
                    rep_options_data = copy.deepcopy(orig_options_data)
                    
                    # Update the options file to dump to a local directory
                    if ('status_files' in rep_options_data['options']):
                        rep_options_data['options']['status_files']['relative_to'] = \
                            'false'
                        last_folder = re.split('/|\\\\',
                                               rep_options_data['options']['status_files']['status_folder'])[-1]
                        rep_options_data['options']['status_files']['status_folder'] = \
                            os.path.join(sim_output_dir,
                                          ('%s_%i_r%i' % (last_folder,
                                                          (k+1), (rep+1))))
                           
                        # Make the status folder if required
                        test_dir = rep_options_data['options']['status_files']['status_folder']
                        
                        if not os.path.isdir(test_dir):
                            os.makedirs(test_dir)
                    
                    # If it is the first force and and the first repeat, update
                    # the options to dump rates
                    if ((k==0) and (rep==0)):
                        rep_options_data['options']['rate_files'] = dict()
                        rep_options_data['options']['rate_files']['relative_to'] = \
                            'false'
                        rep_options_data['options']['rate_files']['file'] = \
                            os.path.join(sim_output_dir,
                                         'rates.json')
                    else:
                        # Clear any existing rate options
                        if ('rate_files' in rep_options_data['options']):
                            del rep_options_data['options']['rate_files']
                    
                    # Create the new options file
                    options_file = os.path.join(sim_input_dir,
                                                ('sim_options_%i_r%i.json' % ((k+1), (rep+1))))
                    
                    with open(options_file, 'w') as f:
                        json.dump(rep_options_data, f, indent=4)
                
                    j = dict()
                    j['relative_to'] = 'False'
                    j['model_file'] = isotonic_model_file
                    j['options_file'] = options_file
                    prot_file_string = os.path.join(sim_input_dir,
                                              ('prot_%i_%i.txt' % ((k+1), (rep+1))))
                    
                    if (fv_struct['pCa'] == 'pCa_50'):
                        test_pCa = c_fit['pCa_50']
                    else:
                        test_pCa = fv_struct['pCa']
                        
                    # Write the protocol
                    df = prot.create_force_control_protocol(
                                            time_step = fv_struct['time_step_s'],
                                            step_pCa = test_pCa,
                                            n_points = int(fv_struct['sim_duration_s'] /
                                                                fv_struct['time_step_s']),
                                            iso_start_s = fv_struct['sim_release_s'],
                                            iso_f = rel_f * isometric_force)
                    
                    if not figures_only:
                        prot.write_protocol_to_file(df, prot_file_string);
                    
                    j['protocol_file'] = prot_file_string

                    # Save the results file
                    j['results_file'] = os.path.join(sim_output_dir,
                                                      ('sim_%i_r%i.txt' % ((k+1), (rep+1))))
        
                    # If required, create an output_handler and add it to
                    # the job
                    if (trace_figures_on == True):
                        print('Need to fix paths for trace figures')
                        exit(1)
                          # Create the structure for the output handler
                        oh = dict()
                        oh['templated_images'] = []
                        tf = dict()
                        tf['relative_to'] = 'this_file'
                        tf['template_file_string'] = os.path.join(
                                                        '..',
                                                        base_dir,
                                                        'template',
                                                        'template_summary.json')
                        tf['output_file_string'] = os.path.join(
                                                        base_dir, fv_struct['sim_folder'],
                                                        'isotonic', 'sim_output',
                                                        ('%i' % (i+1)),
                                                        ('sim_%i_r%i' % ((k+1), (rep+1))))
                        tf['output_image_formats'] = fv_struct['output_image_formats']
                        oh['templated_images'].append(tf)
                    
                        # Now add it to the job, and write it to file
                        j['output_handler_file'] = os.path.join(
                                                    base_dir, fv_struct['sim_folder'],
                                                    'isotonic', 'sim_input',
                                                    ('%i' % (i+1)),
                                                    ('output_handler_sim_%i_r%i.json' %
                                                          ((k+1), (rep+1))))
                    
                        with open(j['output_handler_file'], 'w') as f:
                            json.dump(oh, f, indent=4)      
            
                    isotonic_b['job'].append(j)


    # Now create the batch analysis section
    batch_figs = dict()
    
    # Create the folders for the analysis
    batch_output_dir = str(Path(sim_output_dir).parent)
    
    batch_figs['force_velocity'] = []
    fig = dict()
    fig['relative_to'] = "false"
    fig['results_folder'] = batch_output_dir
    fig['time_release_s'] = fv_struct['sim_release_s']
    fig['fit_time_interval_s'] = fv_struct['fit_time_s']

    if (not 'length_fit_mode' in fv_struct): # fit mode for length traces is not specified, exponential is default
        fig['length_fit_mode'] = 'exponential'
    else:
        fig['length_fit_mode'] = fv_struct['length_fit_mode']


    fig['output_data_file_string'] = os.path.join(batch_output_dir,
                                                  'fv_analysis.xlsx')
    fig['output_image_file'] = os.path.join(batch_output_dir,
                                            'fv_and_power')
    fig['output_image_formats'] = fv_struct['output_image_formats']
    
    if ('formatting' in fv_struct):
        fig['formatting'] = fv_struct['formatting']
    
    batch_figs['force_velocity'].append(fig)

    # Rates
    batch_figs['rates'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = batch_output_dir
    fig['output_image_file'] = os.path.join(batch_output_dir,
                                            'rates')
    fig['output_image_formats'] = fv_struct['output_image_formats']

    if ('formatting' in fv_struct):
        fig['formatting'] = fv_struct['formatting']

    batch_figs['rates'].append(fig)

    # Superposed traces
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = batch_output_dir
    fig['output_image_file'] = os.path.join(batch_output_dir,
                                            'superposed_traces')
    fig['output_image_formats'] = fv_struct['output_image_formats']
    
    if ('formatting' in fv_struct):
        fig['formatting'] = fv_struct['formatting']
        
    batch_figs['superposed_traces'].append(fig)
 
    isotonic_b['batch_figures'] = batch_figs
    
    # Now insert isotonic_b into a full batch structure
    isotonic_batch = dict()
    isotonic_batch['FiberSim_batch'] = isotonic_b

    # Write the batch to file
    base_dir = str(Path(batch_output_dir).parent)
    isotonic_batch_file = os.path.join(base_dir,
                                       'batch_isotonic.json')
    
    with open(isotonic_batch_file, 'w') as f:
        json.dump(isotonic_batch, f, indent=4)
        
    # Now run the isotonic batch
    batch.run_batch(isotonic_batch_file, figures_only=figures_only)

def deduce_freeform_properties(json_analysis_file_string,
                               freeform_struct):
    """ Code runs freeform analysis """
    
    # Potentially switch off simulations
    figures_only = False
    if ('figures_only' in freeform_struct):
        if (freeform_struct['figures_only'] == "True"):
            figures_only = True

    trace_figures_on = True
    if ('trace_figures_on' in freeform_struct):
        if (freeform_struct['trace_figures_on'] == 'False'):
            trace_figures_on = False

    # Load the file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        char_struct = json_data['FiberSim_characterization']
    
    # Pull off the components
    FiberCpp_exe_struct = char_struct['FiberCpp_exe']
    model_struct = char_struct['model']
    
    # Create a batch to run the trials
    freeform_b = dict()
    
    # Turn the FiberCpp_exe into absolute paths because the new instruction
    # file will be in a different place, add to new batch
    cpp_exe = dict()
    if ('relative_to' in FiberCpp_exe_struct):
        if (FiberCpp_exe_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = FiberCpp_exe_struct['relative_to']
        cpp_exe['exe_file'] = str(Path(os.path.join(base_dir,
                                           FiberCpp_exe_struct['exe_file'])).resolve())
    else:
        cpp_exe['exe_file'] = FiberCpp_exe_struct['exe_file']

    freeform_b['FiberCpp_exe'] = cpp_exe

    freeform_b['job'] = []
    
    # Check for half-sarcomere lengths in the pCa_struct
    # If none are specified, create an hsl array from the model file
    if ('hs_lengths' in freeform_struct):
        hs_lengths = freeform_struct['hs_lengths']
    else:
        base_dir = Path(json_analysis_file_string).parent.absolute()
        model_file = os.path.join(base_dir, model_struct['model_files'][0])
        with open(model_file, 'r') as f:
            m = json.load(f)
            hs_lengths = np.array([m['muscle']['initial_hs_length']])
 
    # Set up dir_counter
    dir_counter = 0
    
    # Loop through the model files
    for i, mod_f in enumerate(model_struct['model_files']):

        # Now loop through the half-sarcomere lengths
        for j, hsl in enumerate(hs_lengths):

            # Update dir_counter
            dir_counter = dir_counter + 1
        
            # Create a folder for the sim_input
            if (freeform_struct['relative_to'] == 'this_file'):
                base_dir = Path(json_analysis_file_string).parent.absolute()
            else:
                base_dir = ''
            
            sim_input_dir = os.path.join(base_dir,
                                         freeform_struct['sim_folder'],
                                         'sim_input',
                                         ('%i' % dir_counter))
            
            if not os.path.isdir(sim_input_dir):
                os.makedirs(sim_input_dir)
            
            # Copy the model and options files to the sim_input dir
            # adjusting half-sarcomere lengths as appropriate

            orig_model_file = mod_f
            
            if (model_struct['relative_to'] == 'this_file'):
                base_dir = Path(json_analysis_file_string).parent.absolute()
            else:
                base_dir = ''
            orig_model_file = os.path.join(base_dir, mod_f)
            
            # Adjust hsl by loading model, adjusting hsl and re-writing
            with open(orig_model_file, 'r') as f:
                m = json.load(f)
                m['muscle']['initial_hs_length'] = float(hsl)
                
                # Over-ride m_n if appropriate
                if ('m_n' in freeform_struct):
                    m['thick_structure']['m_n'] = freeform_struct['m_n']
                
            fn = re.split('/|\\\\', orig_model_file)[-1]
            freeform_model_file = os.path.join(sim_input_dir, fn)
            
            with open(freeform_model_file, 'w') as f:
                json.dump(m, f, indent=4)
           
            # Work out the path for the base options file
            orig_options_file = os.path.join(base_dir,
                                                 model_struct['options_file'])

            # Load the options data
            with open(orig_options_file, 'r') as f:
                orig_options_data = json.load(f)

            # Adjust the options struct if we have randomized repeats
            if ('randomized_repeats' in freeform_struct):
                rand_repeats = freeform_struct['randomized_repeats']
                orig_options_data['options']['rand_seed'] = "random"
            else:
                rand_repeats = 1

            # Loop through the protocol_files
            for prot_counter, prot_f in \
                enumerate(freeform_struct['protocol_files']):

                # Loop through the rand_repeats, creating a job for each
                # repeat
                for rep in range(rand_repeats):
                    
                    # Copy the orig_options_struct for local changes
                    # within the rep
                    rep_options_data = copy.deepcopy(orig_options_data)

                    # Update the options file to dump to a local directory
                    if ('status_files' in rep_options_data['options']):
                        rep_options_data['options']['status_files']['status_folder'] = \
                            os.path.join('../../sim_output',
                                         ('%i' % dir_counter),
                                         ('%s_%i_r%i' % (rep_options_data['options']['status_files']['status_folder'],
                                                      (prot_counter + 1), (rep+1))))
                        # Make the status folder if required
                        test_dir = rep_options_data['options']['status_files']['status_folder']
                        if not os.path.isdir(test_dir):
                            os.makedirs(test_dir)

                    if ((prot_counter==0) and (rep == 0)):
                        # If it is the first protocol for the model and the first rep,
                        # update the options file to dump rates
                        rep_options_data['options']['rate_files'] = dict()
                        rep_options_data['options']['rate_files']['relative_to'] = \
                            'this_file'
                        rep_options_data['options']['rate_files']['file'] = \
                            os.path.join('../../sim_output',
                                         ('%i' % dir_counter),
                                         'rates.json')
                    
                    # Create the new options file
                    options_file = os.path.join(
                        Path(sim_input_dir).parent.absolute(),
                        ('%i' % dir_counter),
                        ('sim_options_%i_r%i.json' % (prot_counter+1,
                                                      rep+1)))
        
                    with open(options_file, 'w') as f:
                            json.dump(rep_options_data, f, indent=4)
                    
                    # Copy the protocol file
                    if (freeform_struct['relative_to'] == 'this_file'):
                        base_dir = Path(json_analysis_file_string).parent.absolute()
                    else:
                        base_dir = ''

                    orig_prot_file = os.path.join(base_dir, prot_f)
                    fn = re.split('/|\\\\', orig_prot_file)[-1]
                    freeform_prot_file = os.path.join(sim_input_dir,
                                                      fn)
                    
                    print("opf: %s" % orig_prot_file)
                    print("fpf: %s" % freeform_prot_file)
                    
                    shutil.copyfile(orig_prot_file, freeform_prot_file)

                   # Create the job
                    j = dict()
                    j['relative_to'] = 'False'
                    j['protocol_file'] = freeform_prot_file
                    j['results_file'] = os.path.join(base_dir,
                                                     freeform_struct['sim_folder'],
                                                     'sim_output',
                                                     ('%i' % dir_counter),
                                                     ('sim_prot_%i_r%i.txt' %
                                                      (prot_counter+1, rep+1)))
                    j['model_file'] = freeform_model_file
                    j['options_file'] = options_file
                    
                    # If required, create an output_handler and add it to
                    # the job
                    if (trace_figures_on == True):
                        # Create the structure for the output handler
                        oh = dict()
                        oh['templated_images'] = []
                        tf = dict()
                        tf['relative_to'] = 'this_file'
                        tf['template_file_string'] = os.path.join(
                                                        '..',
                                                        base_dir,
                                                        'template',
                                                        'template_summary.json')
                        tf['output_file_string'] = os.path.join(
                                                        base_dir,
                                                        freeform_struct['sim_folder'],
                                                        'sim_output',
                                                        ('%i' % dir_counter),
                                                        ('sim_prot_%i_r%i' %
                                                         (prot_counter+1, rep+1)))
                        tf['output_image_formats'] = freeform_struct['output_image_formats']
                        oh['templated_images'].append(tf)
                        
                        # Now add it to the job, and write it to file
                        j['output_handler_file'] = os.path.join(
                                                    sim_input_dir,
                                                    ('output_handler_prot_%i_r%i.json' %
                                                        (prot_counter+1, rep+1)))
                        
                        with open(j['output_handler_file'], 'w') as f:
                            json.dump(oh, f, indent=4)        
        
                    freeform_b['job'].append(j)
                
    # Now create the batch analysis section
    batch_figs = dict()
    
    # Rates
    batch_figs['rates'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                         freeform_struct['sim_folder'],
                                         'sim_output')
    fig['output_image_file'] = os.path.join(base_dir,
                                            freeform_struct['sim_folder'],
                                            'sim_output',
                                            'rates')
    fig['output_image_formats'] = freeform_struct['output_image_formats']
    
    if ('formatting' in freeform_struct):
        fig['formatting'] = freeform_struct['formatting']
    
    batch_figs['rates'].append(fig)

    # Superposed traces
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                         freeform_struct['sim_folder'],
                                         'sim_output')
    fig['output_image_file'] = os.path.join(base_dir,
                                            freeform_struct['sim_folder'],
                                            'sim_output',
                                            'superposed_traces')
    fig['output_image_formats'] = freeform_struct['output_image_formats']
    if ('superposed_x_ticks' in freeform_struct):
        fig['superposed_x_ticks'] = freeform_struct['superposed_x_ticks']
    if ('formatting' in freeform_struct):
        fig['formatting'] = freeform_struct['formatting']        
    batch_figs['superposed_traces'].append(fig)
    
    freeform_b['batch_figures'] = batch_figs
                   
    # Now insert isometric_b into a full batch structure
    freeform_batch = dict()
    freeform_batch['FiberSim_batch'] = freeform_b
    
    # Write the batch to file
    base_dir = Path(json_analysis_file_string).parent.absolute()
    freeform_batch_file = os.path.join(base_dir,
                                       freeform_struct['sim_folder'],
                                      'batch_freeform.json')
    with open(freeform_batch_file, 'w') as f:
        json.dump(freeform_batch, f, indent=4)
        
    # Delete existing files from the sim_output folder if running simulations
    if not figures_only:
        sim_output_folder = os.path.join(base_dir,
                                         freeform_struct['sim_folder'],
                                         'sim_output')
        try:
            print('Trying to remove %s' % sim_output_folder)
            shutil.rmtree(sim_output_folder, ignore_errors=True)
        except OSError as e:
            print("Error: %s : %s" % (sim_output_folder, e.strerror))
        
    # Now run the freeform batch
    batch.run_batch(freeform_batch_file, figures_only=figures_only)
    
def characterize_fv_with_pCa_and_isometric_force(json_analysis_file_string,
                                                 fv_characterize_dict):
    """ Runs simulations at given pCa and calculates force-velocity
        properties """
        
    # Potentially switch off simulations
    figures_only = False
    if ('figures_only' in fv_characterize_dict):
        if (fv_characterize_dict['figures_only'] == "True"):
            figures_only = True
        
    # Load the analysis file
    # Load the file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        char_struct = json_data['FiberSim_characterization']
        
    # Create a dict for the force_velocity batch, and fill it in bit by bit
    fv_dict = dict()
        
    # Pull off the exe component
    FiberCpp_exe_struct = char_struct['FiberCpp_exe']
    
    # Turn the FiberCpp_exe into absolute paths because the new instruction
    # file will be in a different place, add to new batch
    cpp_exe = dict()
    cpp_exe['relative_to'] = "False"
    if ('relative_to' in FiberCpp_exe_struct):
        if (FiberCpp_exe_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = FiberCpp_exe_struct['relative_to']
        cpp_exe['exe_file'] = os.path.join(base_dir,
                                           FiberCpp_exe_struct['exe_file'])
    else:
        cpp_exe['exe_file'] = FiberCpp_exe_struct['exe_file']
        
    fv_dict['FiberCpp_exe'] = cpp_exe

    # Turn the model files into absolute paths as well
    model_struct = char_struct['model']
    
    if ('relative_to' in model_struct):
        if (model_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = model_struct['relative_to']
    else:
        base_dir = ''
    
    model_files = []
    for m in model_struct['model_files']:
        model_files.append(os.path.join(base_dir, m))
    
       
    # Check for half-sarcomere lengths in the fv_dict
    # If none are specified, create an hsl array from the model file
    if ('hs_lengths' in fv_dict):
        hs_lengths = fv_dict['hs_lengths']
    else:
        with open(model_files[0], 'r') as f:
            m = json.load(f)
            hs_lengths = np.array([m['muscle']['initial_hs_length']])
             
    # Set up dir_counter
    dir_counter = 0
    
    # Set up your jobs
    fv_dict['job'] = []
    
    # Get the base directory for the simulations
    if ('relative_to' in fv_characterize_dict):
        if (fv_characterize_dict['relative_to'] == 'this_file'):
            base_dir = str(Path(json_analysis_file_string).resolve.absolute().parent)
        else:
            base_dir = fv_characterize_dict['relative_to']
    else:
        base_dir = ''
    base_dir = os.path.join(base_dir, fv_characterize_dict['sim_folder'])
    
    # If you are running simulations, delete the existing structure
    if not figures_only:
        try:
            print('Trying to remove %s' % base_dir)
            shutil.rmtree(base_dir, ignore_errors=True)
        except OSError as e:
            print("Error: %s : %s" % (base_dir, e.strerror))
    
    # Loop through the model files
    for i, mod_f in enumerate(model_struct['model_files']):
        
        # Now loop through the half-sarcomere lengths
        for j, hsl in enumerate(hs_lengths):
            
            # Update dir_counter
            dir_counter = dir_counter + 1
            
            # Create the input and output directories           
            sim_input_dir = os.path.join(base_dir,
                                         'sim_input',
                                         ('%i' % dir_counter))
            if not os.path.isdir(sim_input_dir):
                os.makedirs(sim_input_dir)
                
            sim_output_dir = os.path.join(base_dir,
                                          'sim_output',
                                          ('%i' % dir_counter))
            if not os.path.isdir(sim_output_dir):
                os.makedirs(sim_output_dir)
               
            # Copy the model and options files to the sim_input dir
            # adjusting half-sarcomere length as appropriate
            if (model_struct['relative_to'] == 'this_file'):
                model_dir = Path(json_analysis_file_string).parent.absolute()
                orig_model_file = os.path.join(model_dir, mod_f)
            else:
                orig_model_file = mod_f
            
            # Adjust hsl by loading model, adjusting hsl and re-writing
            with open(orig_model_file, 'r') as f:
                m = json.load(f)
                m['muscle']['initial_hs_length'] = float(hsl)
                
                # Over-ride m_n if appropriate
                if ('m_n' in fv_characterize_dict):
                    m['thick_structure']['m_n'] = fv_characterize_dict['m_n']
                    
            fn = re.split('/|\\\\', orig_model_file)[-1]
            model_file = os.path.join(sim_input_dir, fn)

            with open(model_file, 'w') as f:
                json.dump(m, f, indent=4)
            
            # Work out the path for the base options file
            if (model_struct['relative_to'] == 'this_file'):
                model_dir = Path(json_analysis_file_string).parent.absolute()
                orig_options_file = os.path.join(model_dir,
                                                 model_struct['options_file'])
            else:
                orig_options_file = model_struct['options_file']

            # Load the options data
            with open(orig_options_file, 'r') as f:
                orig_options_data = json.load(f)

            # Adjust the options struct if we have randomized repeats
            if ('randomized_repeats' in fv_characterize_dict):
                rand_repeats = fv_characterize_dict['randomized_repeats']
                orig_options_data['options']['rand_seed'] = "random"
            else:
                rand_repeats = 1
                
            # Loop through the isotonic forces
            for (k, rel_f) in enumerate(fv_characterize_dict['rel_isotonic_forces']):
            
                # Loop through the rand repeats
                for rep in range(rand_repeats):
                    
                    # Copy the orig_options for local changes within the rep
                    rep_options_data = copy.deepcopy(orig_options_data)
                    
                    # Update the options file to dump to a local directory
                    if ('status_files' in rep_options_data['options']):
                        rep_options_data['options']['status_files']['status_folder'] = \
                            os.path.join(sim_output_dir,
                                          ('%s_%i_r%i' % (rep_options_data['options']['status_files']['status_folder'],
                                                          (k+1), (rep+1))))
                    
                    # If it is the first force and and the first repeat, update
                    # the options to dump rates
                    if ((k==0) and (rep==0)):
                        rep_options_data['options']['rate_files'] = dict()
                        rep_options_data['options']['rate_files']['relative_to'] = \
                            'false'
                        rep_options_data['options']['rate_files']['file'] = \
                            os.path.join(sim_output_dir,
                                         'rates.json')
                    else:
                        # Clear any existing rate options
                        if ('rate_files' in rep_options_data['options']):
                            del rep_options_data['options']['rate_files']
                    
                    # Create the new options file
                    options_file = os.path.join(sim_input_dir,
                                                ('sim_options_%i_r%i.json' % ((k+1), (rep+1))))
                    
                    with open(options_file, 'w') as f:
                        json.dump(rep_options_data, f, indent=4)
                
                    j = dict()
                    j['relative_to'] = 'False'
                    j['model_file'] = model_file
                    j['options_file'] = options_file
                    prot_file_string = os.path.join(sim_input_dir,
                                              ('prot_%i_%i.txt' % ((k+1), (rep+1))))
                    
                    print(prot_file_string)
                       
                    # Write the protocol
                    df = prot.create_force_control_protocol(
                                            time_step = fv_characterize_dict['time_step_s'],
                                            step_pCa = fv_characterize_dict['pCa'],
                                            n_points = int(fv_characterize_dict['sim_duration_s'] /
                                                                fv_characterize_dict['time_step_s']),
                                            iso_start_s = fv_characterize_dict['sim_release_s'],
                                            iso_f = rel_f * fv_characterize_dict['isometric_force'])
                    
                    if not figures_only:
                        prot.write_protocol_to_file(df, prot_file_string)
                    
                    j['protocol_file'] = prot_file_string
            
                    # Save the results file
                    j['results_file'] = os.path.join(sim_output_dir,
                                                      ('sim_%i_r%i.txt' % ((k+1), (rep+1))))
            
                    # If required, create an output_handler and add it to
                    # the job
                    # if (trace_figures_on == True):
                    #     print('Need to fix paths for trace figures')
                    #     exit(1)
                    #       # Create the structure for the output handler
                    #     oh = dict()
                    #     oh['templated_images'] = []
                    #     tf = dict()
                    #     tf['relative_to'] = 'this_file'
                    #     tf['template_file_string'] = os.path.join(
                    #                                     '..',
                    #                                     base_dir,
                    #                                     'template',
                    #                                     'template_summary.json')
                    #     tf['output_file_string'] = os.path.join(
                    #                                     base_dir,
                    #                                     fv_characterize_dict['sim_folder'],
                    #                                     'isotonic', 'sim_output',
                    #                                     ('%i' % (i+1)),
                    #                                     ('sim_%i_r%i' % ((k+1), (rep+1))))
                    #     tf['output_image_formats'] = fv_characterize_dict['output_image_formats']
                    #     oh['templated_images'].append(tf)
                    
                    #     # Now add it to the job, and write it to file
                    #     j['output_handler_file'] = os.path.join(
                    #                                 base_dir,
                    #                                 fv_characterize_dict['sim_folder'],
                    #                                 'isotonic', 'sim_input',
                    #                                 ('%i' % (i+1)),
                    #                                 ('output_handler_sim_%i_r%i.json' %
                    #                                       ((k+1), (rep+1))))
                    
                    #     with open(j['output_handler_file'], 'w') as f:
                    #         json.dump(oh, f, indent=4)      
            
                    fv_dict['job'].append(j)          
    
    # Now create the batch analysis section
    batch_figs = dict()
    
    # Create the folders for the analysis
    batch_output_dir = str(Path(sim_output_dir).parent)
    
    batch_figs['force_velocity'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = batch_output_dir
    fig['time_release_s'] = fv_characterize_dict['sim_release_s']
    fig['fit_time_interval_s'] = fv_characterize_dict['fit_time_s']

    if (not 'length_fit_mode' in fv_characterize_dict): # fit mode for length traces is not specified, exponential is default
        fig['length_fit_mode'] = 'exponential'
    else:
        fig['length_fit_mode'] = fv_characterize_dict['length_fit_mode']


    fig['output_data_file_string'] = os.path.join(batch_output_dir,
                                                  'fv_analysis.xlsx')
    fig['output_image_file'] = os.path.join(batch_output_dir,
                                            'fv_and_power')
    fig['output_image_formats'] = fv_characterize_dict['output_image_formats']
    
    if ('formatting' in fv_characterize_dict):
        fig['formatting'] = fv_characterize_dict['formatting']
    
    batch_figs['force_velocity'].append(fig)

    # Rates
    batch_figs['rates'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = batch_output_dir
    fig['output_image_file'] = os.path.join(batch_output_dir,
                                            'rates')
    fig['output_image_formats'] = fv_characterize_dict['output_image_formats']

    if ('formatting' in fv_characterize_dict):
        fig['formatting'] = fv_characterize_dict['formatting']

    batch_figs['rates'].append(fig)

    # Superposed traces
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = batch_output_dir
    fig['output_image_file'] = os.path.join(batch_output_dir,
                                            'superposed_traces')
    fig['output_image_formats'] = fv_characterize_dict['output_image_formats']
    
    if ('formatting' in fv_characterize_dict):
        fig['formatting'] = fv_characterize_dict['formatting']
        
    batch_figs['superposed_traces'].append(fig)
 
    fv_dict['batch_figures'] = batch_figs
    
    # Now insert isotonic_b into a full batch structure
    isotonic_batch = dict()
    isotonic_batch['FiberSim_batch'] = fv_dict

    # Write the batch to file
    base_dir = str(Path(batch_output_dir).parent)
    isotonic_batch_file = os.path.join(base_dir,
                                       'batch_isotonic.json')
    
    with open(isotonic_batch_file, 'w') as f:
        json.dump(isotonic_batch, f, indent=4)
        
    # Now run the isotonic batch
    batch.run_batch(isotonic_batch_file, figures_only=figures_only)
        