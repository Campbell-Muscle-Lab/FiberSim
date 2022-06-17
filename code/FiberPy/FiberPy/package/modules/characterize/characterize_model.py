# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 12:47:26 2022

@author: kscamp3
"""
import os
import json
import shutil

import numpy as np
import pandas as pd

from pathlib import Path

from ..protocols import protocols as prot
from ..batch import batch

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
        anal_struct = json_data['FiberSim_analysis']
   
    # Pull of the analysis tasks
    anal_struct = anal_struct['analysis']
    for an in anal_struct:
        if (an['type'] == 'pCa_length_control'):
            deduce_pCa_length_control_properties(json_analysis_file_string,
                                                 pCa_struct = an)
        if (an['type'] == 'force_velocity'):
            deduce_fv_properties(json_analysis_file_string,
                                 fv_struct = an)
        
        if (an['type'] == 'freeform'):
            deduce_freeform_properties(json_analysis_file_string,
                                       freeform_struct = an)
            
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
        anal_struct = json_data['FiberSim_analysis']
    
    # Pull off the components
    FiberCpp_exe_struct = anal_struct['FiberCpp_exe']
    model_struct = anal_struct['model']
    
    # Create an isometric batch to run the isomtric test
    pCa_lc_b = dict()
    
    # Turn the FiberCpp_exe into absolute paths because the new instruction
    # file will be in a different place
    if (FiberCpp_exe_struct['relative_to'] == 'this_file'):
        base_dir = Path(json_analysis_file_string).parent.absolute()
        FiberCpp_exe_struct['relative_to'] = 'False'
        FiberCpp_exe_struct['exe_file'] = \
            os.path.join(base_dir, FiberCpp_exe_struct['exe_file'])
        pCa_lc_b['FiberCpp_exe'] = FiberCpp_exe_struct

    pCa_lc_b['job'] = []
    
    # Check for half-sarcomere lengths in the pCa_struct
    # If none are specified, create an hsl array from the model file
    if ('hs_lengths' in pCa_struct):
        hs_lengths = pCa_struct['hs_lengths']
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
            if (pCa_struct['relative_to'] == 'this_file'):
                base_dir = Path(json_analysis_file_string).parent.absolute()
                sim_input_dir = os.path.join(base_dir,
                                             pCa_struct['sim_folder'],
                                             'sim_input',
                                             ('%i' % dir_counter))
                if not os.path.isdir(sim_input_dir):
                    os.makedirs(sim_input_dir)
               
            # Copy the model and options files to the sim_input dir
            # adjusting half-sarcomere length as appropriate
            if (model_struct['relative_to'] == 'this_file'):
                base_dir = Path(json_analysis_file_string).parent.absolute()
                
                orig_model_file = os.path.join(base_dir, mod_f)
                
                # Adjust hsl by loading model, adjusting hsl and re-writing
                with open(orig_model_file, 'r') as f:
                    m = json.load(f)
                    m['muscle']['initial_hs_length'] = float(hsl)
                fn = orig_model_file.split('/')[-1]
                iso_model_file = os.path.join(sim_input_dir, fn)

                with open(iso_model_file, 'w') as f:
                    json.dump(m, f, indent=4)
                
                # Now copy the options file
                orig_options_file = os.path.join(base_dir, model_struct['options_file'])
                fn = orig_options_file.split(os.sep)[-1]
                iso_options_file = os.path.join(sim_input_dir, fn)
                shutil.copyfile(orig_options_file, iso_options_file)
           
            # Loop through the pCa values
            for pCa_counter,pCa in enumerate(pCa_struct['pCa_values']):
                
                if (pCa_counter==0):
                    # Update the options file to dump rates
                    with open(iso_options_file, 'r') as f:
                        json_data = json.load(f)
                        json_data['options']['rate_files'] = dict()
                        json_data['options']['rate_files']['relative_to'] = 'this_file'
                        json_data['options']['rate_files']['file'] = \
                            os.path.join('../../sim_output',
                                         ('%i' % dir_counter),
                                         'rates.txt')
                    iso_options_file_rates = os.path.join(
                        Path(iso_options_file).parent.absolute(),
                        'sim_options_rates.json')
    
                    with open(iso_options_file_rates, 'w') as f:
                        json.dump(json_data, f, indent=4)
    
                # Set the delta_hsl vector
                n_points = int(pCa_struct['sim_duration_s'] /
                               pCa_struct['time_step_s'])
                delta_hsl = np.zeros(n_points)
                mode_vector = -2 * np.ones(n_points)
                
                # Add in k_tr if required
                if ('k_tr_start_s' in pCa_struct):
                    # Calculate some stuff for the k_tr
                    k_tr_start_ind = int(pCa_struct['k_tr_start_s'] /
                                         pCa_struct['time_step_s'])
                    k_tr_stop_ind = int((pCa_struct['k_tr_start_s'] + 
                                             pCa_struct['k_tr_duration_s']) /
                                        pCa_struct['time_step_s'])
                    k_tr_ramp_points = int(pCa_struct['k_tr_ramp_s'] /
                                           pCa_struct['time_step_s'])
                    ramp_inc = pCa_struct['k_tr_magnitude_nm'] / \
                                    float(k_tr_ramp_points)
                    # Set the k_tr_shortening
                    vi = np.arange(k_tr_start_ind,
                                   k_tr_start_ind + k_tr_ramp_points, 1)
                    delta_hsl[vi] = -ramp_inc
                    # Set the k_tr_re-stretch
                    vi = np.arange(k_tr_stop_ind,
                                   k_tr_stop_ind + k_tr_ramp_points, 1)
                    delta_hsl[vi] = ramp_inc
                    # Set the mode
                    vi = np.arange(k_tr_start_ind,
                                   k_tr_stop_ind + k_tr_ramp_points, 1)
                    mode_vector[vi] = -1
                
                # Add in user-defined delta_hsl
                if ('user_defined_dhsl' in pCa_struct):
                    dhsl_file = os.path.join(base_dir,
                                             pCa_struct['user_defined_dhsl'])
                    dhsl = pd.read_csv(dhsl_file)
                    delta_hsl = dhsl['dhsl'].to_numpy()
                    
                # Create a length control protocol and write to file
                df = prot.create_length_control_protocol(
                        time_step = pCa_struct['time_step_s'],
                        step_pCa = pCa,
                        n_points = n_points,
                        delta_hsl = delta_hsl,
                        mode_vector = mode_vector)
                prot_file_string = os.path.join(sim_input_dir,
                                                'prot_pCa_%.0f.txt' % (10*pCa))
                prot.write_protocol_to_file(df, prot_file_string)
                
                # Create the job
                j = dict()
                j['relative_to'] = 'False'
                j['protocol_file'] = prot_file_string
                j['results_file'] = os.path.join(base_dir,
                                                 pCa_struct['sim_folder'],
                                                 'sim_output',
                                                 ('%i' % dir_counter),
                                                 'sim_pCa_%.0f.txt' % (10*pCa))
                j['model_file'] = iso_model_file
                if (pCa_counter == 0):
                    j['options_file'] = iso_options_file_rates
                else:
                    j['options_file'] = iso_options_file
                
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
                                                    pCa_struct['sim_folder'],
                                                    'sim_output',
                                                    ('%i' % dir_counter),
                                                    'sim_pCa_%.0f' % (10*pCa))
                    tf['output_image_formats'] = pCa_struct['output_image_formats']
                    oh['templated_images'].append(tf)
                    
                    # Now add it to the job, and write it to file
                    j['output_handler_file'] = os.path.join(
                                                sim_input_dir,
                                                'output_handler_pCa_%.0f.json' %
                                                    (10*pCa))
                    
                    with open(j['output_handler_file'], 'w') as f:
                        json.dump(oh, f, indent=4)        
    
                pCa_lc_b['job'].append(j)
    
    # Now create the analysis section
    batch_figs = dict()

    # pCa curves
    batch_figs['pCa_curves'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                          pCa_struct['sim_folder'],
                                          'sim_output')
    fig['data_field'] = 'force'
    fig['output_data_file_string'] = os.path.join(base_dir,
                                                  pCa_struct['sim_folder'],
                                                  'sim_output',
                                                  'pCa_analysis.xlsx')
    fig['output_image_file'] = os.path.join(base_dir,
                                                   pCa_struct['sim_folder'],
                                                   'sim_output',
                                                   'force_pCa')
    fig['output_image_formats'] = pCa_struct['output_image_formats']
    if ('formatting' in pCa_struct):
        fig['formatting'] = pCa_struct['formatting']
    else:
        fig['formatting'] = dict()
        fig['formatting']['y_axis_label'] = 'Force (N m$^{\\mathregular{-2}}$)'
    batch_figs['pCa_curves'].append(fig)

    # pCa curves - NORMALIZED
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                          pCa_struct['sim_folder'],
                                          'sim_output')
    fig['data_field'] = 'force'
    fig['output_image_file'] = os.path.join(base_dir,
                                                   pCa_struct['sim_folder'],
                                                   'sim_output',
                                                   'force_pCa_normalized')
    fig['output_image_formats'] = pCa_struct['output_image_formats']
    fig['formatting'] = dict()

    fig['formatting']['y_axis_label'] = 'Normalized \n force'
    fig['formatting']['y_normalized_to_max'] = 'True'
    fig['formatting']['y_label_pad'] = 20

    batch_figs['pCa_curves'].append(fig)

    # Rates
    batch_figs['rates'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                         pCa_struct['sim_folder'],
                                         'sim_output')
    fig['output_image_file'] = os.path.join(base_dir,
                                            pCa_struct['sim_folder'],
                                            'sim_output',
                                            'rates')
    fig['output_image_formats'] = pCa_struct['output_image_formats']
    
    if ('formatting' in pCa_struct):
        fig['formatting'] = pCa_struct['formatting']
    
    batch_figs['rates'].append(fig)

    # Superposed traces
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                         pCa_struct['sim_folder'],
                                         'sim_output')
    fig['output_image_file'] = os.path.join(base_dir,
                                            pCa_struct['sim_folder'],
                                            'sim_output',
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
        fig['results_folder'] = os.path.join(base_dir,
                                         pCa_struct['sim_folder'],
                                         'sim_output')
        fig['output_data_file_string'] = os.path.join(
                                             base_dir,
                                             pCa_struct['sim_folder'],
                                             'sim_output',
                                             'k_tr_analysis.xlsx')
        fig['output_image_file'] = os.path.join(base_dir,
                                            pCa_struct['sim_folder'],
                                            'sim_output',
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
    base_dir = Path(json_analysis_file_string).parent.absolute()
    pCa_lc_batch_file = os.path.join(base_dir,
                                     pCa_struct['sim_folder'],
                                     'batch_isometric.json')
    with open(pCa_lc_batch_file, 'w') as f:
        json.dump(pCa_lc_batch, f, indent=4)
        
    # Now run the isometric batch
    batch.run_batch(pCa_lc_batch_file, figures_only=figures_only)
    

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
        anal_struct = json_data['FiberSim_analysis']
    
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

    isometric_b['job'] = []

    # Loop through the model files
    for i, mod_f in enumerate(model_struct['model_files']):

        # Create a folder for the sim_input
        if (fv_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
            sim_input_dir = os.path.join(base_dir, fv_struct['sim_folder'],
                                         'isometric', 'sim_input',
                                         ('%i' % (i+1)))
            if not os.path.isdir(sim_input_dir):
                os.makedirs(sim_input_dir)
            
        # Copy the model and options files to the sim_input dir
        if (model_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
            
            orig_model_file = os.path.join(base_dir, mod_f)
            fn = orig_model_file.split('/')[-1]
            iso_model_file = os.path.join(sim_input_dir, fn)
            shutil.copyfile(orig_model_file, iso_model_file)
            
            orig_options_file = os.path.join(base_dir, model_struct['options_file'])
            fn = orig_options_file.split(os.sep)[-1]
            iso_options_file = os.path.join(sim_input_dir, fn)
            shutil.copyfile(orig_options_file, iso_options_file)
            
        
            # Create a length control protocol and write it to file
            df = prot.create_length_control_protocol(
                                    time_step = fv_struct['time_step_s'],
                                    step_pCa = fv_struct['pCa'],
                                    n_points = int(fv_struct['sim_duration_s'] /
                                                        fv_struct['time_step_s']))
            prot_file_string = os.path.join(sim_input_dir,
                                            'prot_iso_pCa_%.0f.txt' % (10 * fv_struct['pCa']))
            prot.write_protocol_to_file(df, prot_file_string);
 
        j = dict()
        j['relative_to'] = 'False'
        j['protocol_file'] = prot_file_string
        j['results_file'] = os.path.join(base_dir, fv_struct['sim_folder'],
                                         'isometric', 'sim_output',
                                         ('%i' % (i+1)),
                                         'sim_pCa_%.0f.txt' % (10 * fv_struct['pCa']))
        j['model_file'] = iso_model_file
        j['options_file'] = iso_options_file

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
                                            base_dir, fv_struct['sim_folder'],
                                            'isometric', 'sim_output',
                                            ('%i' % (i+1)),
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
    
    # Write the batch to file
    base_dir = Path(json_analysis_file_string).parent.absolute()
    isometric_batch_file = os.path.join(base_dir, fv_struct['sim_folder'],
                                        'isometric',
                                        'batch_isometric.json')
    with open(isometric_batch_file, 'w') as f:
        json.dump(isometric_batch, f, indent=4)
        
    # Now run the isometric batch
    batch.run_batch(isometric_batch_file, figures_only=figures_only)

    # Switch to isotonic stuff
    # First create the isotonic batch dict
    isotonic_b = dict()
    isotonic_b['FiberCpp_exe'] = FiberCpp_exe_struct
    isotonic_b['job'] = []
    
    # Now cycle thought the isometric jobs, generating an isotonic suite
    # for each one
    for (i,r) in enumerate(isometric_batch['FiberSim_batch']['job']):

        # Pull off the isometric force
        sim_data = pd.read_csv(r['results_file'], sep='\t')
        isometric_force = sim_data['force'].iloc[-50:-1].mean() # take the mean force over last 50 points

        # Deduce the sim input dir
        base_dir = Path(json_analysis_file_string).parent.absolute()
        sim_input_dir = os.path.join(base_dir, fv_struct['sim_folder'],
                                     'isotonic', 'sim_input',
                                     ('%i' % (i+1)))
        # Make it if required
        if not os.path.isdir(sim_input_dir):
            os.makedirs(sim_input_dir)
        
        orig_model_file = r['model_file']
        fn = orig_model_file.split('\\')[-1]
        isotonic_model_file = os.path.join(sim_input_dir, fn)
        shutil.copyfile(orig_model_file, isotonic_model_file)
        
        orig_options_file = r['options_file']
        fn = orig_options_file.split('\\')[-1]
        isotonic_options_file = os.path.join(sim_input_dir, fn)
        shutil.copyfile(orig_options_file, isotonic_options_file)

        # Cycle through the isotonic forces

        for (k, rel_f) in enumerate(fv_struct['rel_isotonic_forces']):

            if (k==0):
                # Update the options file to dump rates
                with open(isotonic_options_file, 'r') as f:
                    json_data = json.load(f)
                    json_data['options']['rate_files'] = dict()
                    json_data['options']['rate_files']['relative_to'] = 'this_file'
                    json_data['options']['rate_files']['file'] = \
                        os.path.join('../../sim_output',
                                     ('%i' % (i+1)),
                                     'rates.txt')
                isotonic_options_file_rates = os.path.join(
                    Path(isotonic_options_file).parent.absolute(),
                    'sim_options_rates.json')

                with open(isotonic_options_file_rates, 'w') as f:
                    json.dump(json_data, f, indent=4)
            
            j = dict()
            j['relative_to'] = 'False'
            j['model_file'] = isotonic_model_file
            if (k==0):
                j['options_file'] = isotonic_options_file_rates
            else:
                j['options_file'] = isotonic_options_file
            prot_file_string = os.path.join(sim_input_dir,
                                      ('prot_%i.txt' % (k+1)))
            df = prot.create_force_control_protocol(
                                    time_step = fv_struct['time_step_s'],
                                    step_pCa = fv_struct['pCa'],
                                    n_points = int(fv_struct['sim_duration_s'] /
                                                        fv_struct['time_step_s']),
                                    iso_start_s = fv_struct['sim_release_s'],
                                    iso_f = rel_f * isometric_force)
            prot.write_protocol_to_file(df, prot_file_string);
            j['protocol_file'] = prot_file_string
            j['results_file'] = os.path.join(base_dir, fv_struct['sim_folder'],
                                              'isotonic', 'sim_output',
                                              ('%i' % (i+1)),
                                              ('sim_%i.txt' % (k+1)))

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
                                                base_dir, fv_struct['sim_folder'],
                                                'isotonic', 'sim_output',
                                                ('%i' % (i+1)),
                                                ('sim_%i' % (k+1)))
                tf['output_image_formats'] = fv_struct['output_image_formats']
                oh['templated_images'].append(tf)
            
                # Now add it to the job, and write it to file
                j['output_handler_file'] = os.path.join(
                                            base_dir, fv_struct['sim_folder'],
                                            'isotonic', 'sim_input',
                                            ('%i' % (i+1)),
                                            ('output_handler_sim_%i.json' % (k+1)))
            
                with open(j['output_handler_file'], 'w') as f:
                    json.dump(oh, f, indent=4)      

            isotonic_b['job'].append(j)

    # Now create the batch analysis section
    batch_figs = dict()
    batch_figs['force_velocity'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                          fv_struct['sim_folder'],
                                          'isotonic',
                                          'sim_output')
    fig['time_release_s'] = fv_struct['sim_release_s']
    fig['fit_time_interval_s'] = fv_struct['fit_time_s']

    if (not 'length_fit_mode' in fv_struct): # fit mode for length traces is not specified, exponential is default
        fig['length_fit_mode'] = 'exponential'
    else:
        fig['length_fit_mode'] = fv_struct['length_fit_mode']


    fig['output_data_file_string'] = os.path.join(base_dir,
                                                  fv_struct['sim_folder'],
                                                  'isotonic',
                                                  'sim_output',
                                                  'fv_analysis.xlsx')
    fig['output_image_file'] = os.path.join(base_dir,
                                            fv_struct['sim_folder'],
                                            'isotonic',
                                            'sim_output',
                                            'fv_and_power')
    fig['output_image_formats'] = fv_struct['output_image_formats']
    
    if ('formatting' in fv_struct):
        fig['formatting'] = fv_struct['formatting']
    
    batch_figs['force_velocity'].append(fig)

    # Rates
    batch_figs['rates'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                         fv_struct['sim_folder'],
                                         'isotonic',
                                         'sim_output')
    fig['output_image_file'] = os.path.join(base_dir,
                                            fv_struct['sim_folder'],
                                            'isotonic',
                                            'sim_output',
                                            'rates')
    fig['output_image_formats'] = fv_struct['output_image_formats']

    if ('formatting' in fv_struct):
        fig['formatting'] = fv_struct['formatting']

    batch_figs['rates'].append(fig)

    # Superposed traces
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(base_dir,
                                         fv_struct['sim_folder'],
                                         'isotonic',
                                         'sim_output')
    fig['output_image_file'] = os.path.join(base_dir,
                                            fv_struct['sim_folder'],
                                            'isotonic',
                                            'sim_output',
                                            'superposed_traces')
    fig['output_image_formats'] = fv_struct['output_image_formats']
    
    if ('formatting' in fv_struct):
        fig['formatting'] = fv_struct['formatting']
        
    batch_figs['superposed_traces'].append(fig)
 
    isotonic_b['batch_figures'] = batch_figs
    
    # Now insert isotonic_b into a full batch structure
    isotonic_batch = dict()
    isotonic_batch['FiberSim_batch'] = isotonic_b

    base_dir = Path(json_analysis_file_string).parent.absolute()
    isotonic_batch_file = os.path.join(base_dir, fv_struct['sim_folder'],
                                        'isotonic',
                                        'batch_isotonic.json')
    
    with open(isotonic_batch_file, 'w') as f:
        json.dump(isotonic_batch, f, indent=4)
        
    # Now run the isotonic batch
    batch.run_batch(isotonic_batch_file, figures_only=figures_only)

def deduce_freeform_properties(json_analysis_file_string,
                               freeform_struct):
    """ Code runs freeform analysis """
    
    print(freeform_struct)
    
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
        anal_struct = json_data['FiberSim_analysis']
    
    # Pull off the components
    FiberCpp_exe_struct = anal_struct['FiberCpp_exe']
    model_struct = anal_struct['model']
    
    # Create a batch to run the trials
    freeform_b = dict()
    
    # Turn the FiberCpp_exe into absolute paths because the new instruction
    # file will be in a different place
    if (FiberCpp_exe_struct['relative_to'] == 'this_file'):
        base_dir = Path(json_analysis_file_string).parent.absolute()
        FiberCpp_exe_struct['relative_to'] = 'False'
        FiberCpp_exe_struct['exe_file'] = \
            os.path.join(base_dir, FiberCpp_exe_struct['exe_file'])
        freeform_b['FiberCpp_exe'] = FiberCpp_exe_struct

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
                sim_input_dir = os.path.join(base_dir,
                                             freeform_struct['sim_folder'],
                                             'sim_input',
                                             ('%i' % dir_counter))
                if not os.path.isdir(sim_input_dir):
                    os.makedirs(sim_input_dir)
                
            # Copy the model and options files to the sim_input dir
            # adjusting half-sarcomere lengths as appropriate
       
            if (model_struct['relative_to'] == 'this_file'):
                base_dir = Path(json_analysis_file_string).parent.absolute()
                
                orig_model_file = os.path.join(base_dir, mod_f)
                
                # Adjust hsl by loading model, adjusting hsl and re-writing
                with open(orig_model_file, 'r') as f:
                    m = json.load(f)
                    m['muscle']['initial_hs_length'] = float(hsl)
                    
                    # Over-ride m_n if appropriate
                    if ('m_n' in freeform_struct):
                        m['thick_structure']['m_n'] = freeform_struct['m_n']
                    
                fn = orig_model_file.split('/')[-1]
                freeform_model_file = os.path.join(sim_input_dir, fn)
                
                with open(freeform_model_file, 'w') as f:
                    json.dump(m, f, indent=4)
                
                # Now copy the options file
                orig_options_file = os.path.join(base_dir,
                                                 model_struct['options_file'])
                fn = orig_options_file.split(os.path.sep)[-1]
                freeform_options_file = os.path.join(sim_input_dir, fn)
                
                print(fn)
                print(freeform_options_file)
                
                shutil.copyfile(orig_options_file, freeform_options_file)
           
            # Loop through the protocol_files
            for prot_counter, prot_f in \
                enumerate(freeform_struct['protocol_files']):
                
                if (prot_counter==0):
                    # Update the options file to dump rates
                    with open(freeform_options_file, 'r') as f:
                        json_data = json.load(f)
                        json_data['options']['rate_files'] = dict()
                        json_data['options']['rate_files']['relative_to'] = \
                            'this_file'
                        json_data['options']['rate_files']['file'] = \
                            os.path.join('../../sim_output',
                                         ('%i' % dir_counter),
                                         'rates.txt')
                    freeform_options_file_rates = os.path.join(
                        Path(freeform_options_file).parent.absolute(),
                        'sim_options_rates.json')
    
                    with open(freeform_options_file_rates, 'w') as f:
                        json.dump(json_data, f, indent=4)
                
                # Copy the protocol file
                orig_prot_file = os.path.join(base_dir, prot_f)
                fn = orig_prot_file.split('/')[-1]
                freeform_prot_file = os.path.join(sim_input_dir, fn)
                shutil.copyfile(orig_prot_file, freeform_prot_file)
                
                # Create the job
                j = dict()
                j['relative_to'] = 'False'
                j['protocol_file'] = freeform_prot_file
                j['results_file'] = os.path.join(base_dir,
                                                 freeform_struct['sim_folder'],
                                                 'sim_output',
                                                 ('%i' % dir_counter),
                                                 ('sim_prot_%i.txt' % prot_counter))
                j['model_file'] = freeform_model_file
                if (prot_counter == 0):
                    j['options_file'] = freeform_options_file_rates
                else:
                    j['options_file'] = freeform_options_file
                
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
                                                    ('sim_prot_%i' % prot_counter))
                    tf['output_image_formats'] = freeform_struct['output_image_formats']
                    oh['templated_images'].append(tf)
                    
                    # Now add it to the job, and write it to file
                    j['output_handler_file'] = os.path.join(
                                                sim_input_dir,
                                                ('output_handler_prot_%i.json' %
                                                    prot_counter))
                    
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
        
    # Now run the freeform batch
    batch.run_batch(freeform_batch_file, figures_only=figures_only)