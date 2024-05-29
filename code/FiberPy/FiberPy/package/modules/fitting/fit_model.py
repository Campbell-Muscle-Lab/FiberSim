# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:36:17 2024

@author: Campbell
"""

import os
import json
import shutil
import copy
import subprocess

import numpy as np
import pandas as pd

from scipy.optimize import minimize

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from ..batch import batch

from ..characterize import characterize_model

def fit_model(json_analysis_file_string):
    """ Code takes a setup file with a model/fitting section, and
        tries to fit the model to data """
    
    # Check the analysis file
    if (not json_analysis_file_string):
        print('fitting: fit_model: no analysis file specified')
        exit(1)
        
    # Load it
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        model_struct = json_data['FiberSim_setup']['model']
        fitting_struct = model_struct['fitting']

    # Deduce the number of parameters
    
    # In the simplest case, the number of parameters is defined
    # by the number of 'simple' adjustments
    no_of_parameters = len(fitting_struct['adjustments'])
    
    # We also need to check for base_variants
    if ('base_variants' in fitting_struct):
        base_variants = fitting_struct['base_variants']
        for (i, bv) in enumerate(base_variants):
            adjustments = bv['adjustments']
            for (j, adj) in enumerate(adjustments):
                if not ('multipliers' in adj):
                    no_of_parameters = no_of_parameters + 1

    # Now set p    
    if ('initial_guess' in fitting_struct):
        p = np.asarray(fitting_struct['initial_guess'])
    else:
        p = 0.5 * np.ones(no_of_parameters)
        
    # Set up for fitting
    
    # Set the progress folder
    if (model_struct['relative_to'] == 'this_file'):
        model_base_dir = Path(json_analysis_file_string).parent
    else:
        model_base_dir = model_struct['relative_to']
    progress_dir = os.path.join(model_base_dir,
                                   fitting_struct['progress_folder'])
    progress_dir = str(Path(progress_dir).resolve().absolute())
    
    # Clean progress_dir
    try:
        print('Trying to clean: %s' % progress_dir)
        shutil.rmtree(progress_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s' % (progress_dir, e.strerror))
        
    # Check the working dir is there
    if not os.path.isdir(progress_dir):
        os.makedirs(progress_dir)

    # And now the current folder
    current_dir = os.path.join(model_base_dir, fitting_struct['current_folder'])
    current_dir = str(Path(current_dir).resolve().absolute())

    # Clean current dir
    try:
        print('Trying to clean: %s', current_dir)
        shutil.rmtree(current_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s' % (current_dir, e.strerror))

    # Check the directory isthere
    if not os.path.isdir(current_dir):
        os.makedirs(current_dir)
    
    # Make the dictionary
    progress_data = dict()
    progress_data['iteration'] = 1
    progress_data['lowest_error'] = np.inf
    progress_data['best_p_vector'] = p
    progress_data['error_values'] = []
    progress_data['progress_folder'] = progress_dir
    progress_data['progress_file_string'] = os.path.join(progress_dir,
                                                         'progress.xlsx')
    progress_data['best_file_string'] = os.path.join(progress_dir,
                                                     'best.xlsx')
    progress_data['current_folder'] = current_dir
        
    if (fitting_struct['single_run'] == 'True'):
        print('Single run')
        worker(p, json_analysis_file_string, progress_data)
        print('Finished single run')
        return
    
    
    # Set up the optimizer
    if (fitting_struct['optimizer'] == 'particle_swarm'):
        
        pso(worker, p, json_analysis_file_string, progress_data)
   
    else:
        
        # Create the bounds
        bnds = []
        for i in range(no_of_parameters):
            bnds.append(tuple([0, 1]))
        bnds = tuple(bnds)
        
        res = minimize(worker, p,
                       args=(json_analysis_file_string, progress_data),
                       method=fitting_struct['optimizer'],
                       bounds=bnds)
        
def pso(worker, p_vector, json_analysis_file_string, progress_data,
        f_particles = 2, bounds = [0, 1],
        inertia = 0.5, w_self = 0.5, w_family = 0.5,
        initial_vel = 0.0,
        vel_bounds = [0.0, 0.2],
        jitter = 0.03,
        throw = 20):
    
    # Set the number of particles
    n_particles = round(f_particles * len(p_vector))
    
    # Set the initial values
    rng = np.random.default_rng()
    x = rng.random([n_particles, len(p_vector)])
    x[0,:] = p_vector

    # Set the initial values
    global_best_value = np.Inf
    global_best_x = np.NaN * np.ones(len(p_vector))
    particle_best_value = np.Inf * np.ones(n_particles)
    particle_best_x = np.NaN * np.ones([n_particles, len(p_vector)])
    v = initial_vel * np.ones([n_particles, len(p_vector)])
    
    for iter in range(100):
        for i in range(n_particles):

            particle_value = worker(x[i,:],
                                      json_analysis_file_string,
                                      progress_data)
            
            if (particle_value < particle_best_value[i]):
                particle_best_value[i] = particle_value
                particle_best_x[i,:] = copy.deepcopy(x[i,:])
            
            if (particle_best_value[i] < global_best_value):
                global_best_value = particle_best_value[i]
                global_best_x = copy.deepcopy(x[i,:])
                
        for i in range(n_particles):
            print("particle_best_value[%i]: %g" % (i, particle_best_value[i]))
            print(particle_best_x[i,:])
                
        # Update
        for i in range(n_particles):
            for j in range(len(p_vector)):
                if (x[i,j] == bounds[0]):
                    v[i,j] = 0
                if (x[i,j] == bounds[-1]):
                    v[i,j] = 0
                
                v[i,j] = inertia*v[i,j] + \
                    w_self * rng.random() * (particle_best_x[i,j] - x[i,j]) + \
                    w_family * rng.random() * (global_best_x[j] - x[i,j])
                
                if (v[i,j] > vel_bounds[-1]):
                    v[i,j] = vel_bounds[-1]
                if (v[i,j] < -vel_bounds[-1]):
                    v[i,j]  = -vel_bounds[-1]
                
                x[i,j] = x[i,j] + v[i,j]
                
                # Add in some jitter
                x[i,j] = x[i,j] + jitter * (rng.random() - 0.5)
                
                if (x[i,j] < bounds[0]):
                    x[i,j] = bounds[0]
                if (x[i,j] > bounds[-1]):
                    x[i,j] = bounds[-1]
                    
        # Throw
        if ((iter % throw) == 0):
            # Throw out 1/3 of particles
            for i in range(round(n_particles / 3)):
                for j in range(len(p_vector)):
                    x[i,j] = rng.random()
                    v[i,j] = 0
                particle_best_value[i] = np.Inf
                particle_best_x[i,:] = x[i,:]
    
    print('Done')

    
def worker(p_vector, json_analysis_file_string, progress_data):
    """ Code launches a simulation with parameter multipliers set
        by p_vector and returns a single error value """
        
    
    # Open the analysis file
    with open(json_analysis_file_string, 'r') as f:
        orig_setup = json.load(f)
        
    # Copy the dict
    new_setup = copy.deepcopy(orig_setup)
    
    # Pull out the model_struct
    model_struct = orig_setup['FiberSim_setup']['model']
    
    # Deduce the base directory
    if (model_struct['relative_to'] == 'this_file'):
        model_base_dir = Path(json_analysis_file_string).parent.absolute()
    else:
        model_base_dir = model_struct['relative_to']

    # Set the working directory
    working_dir = os.path.join(model_base_dir,
                               model_struct['fitting']['working_folder'])
    
    # Clean the working directory
    working_dir = os.path.join(model_base_dir,
                               model_struct['fitting']['working_folder'])
    try:
        print('Trying to clean: %s' % working_dir)
        shutil.rmtree(working_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s' % (working_dir, e.strerror))
        
    # Check the working dir is there
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
        
    # Move the model file and the simulations to the working dir
    file_names = [model_struct['fitting']['base_model'],
                  model_struct['options_file']]
    for fn in file_names:
        orig_fn = os.path.join(model_base_dir, fn)
        with open(orig_fn, 'r') as f:
            orig_data = json.load(f)
        new_fn = os.path.join(working_dir, fn)
        with open(new_fn, 'w') as f:
            json.dump(orig_data, f, indent=4)
        
    # Rename the fitting key
    new_setup['FiberSim_setup']['model']['manipulations'] = \
        new_setup['FiberSim_setup']['model'].pop('fitting')
        
    # Remove keys we don't need
    for k in ['working_folder', 'progress_folder','Python_objective_call']:
        new_setup['FiberSim_setup']['model']['manipulations'].pop(k)
        
    # Set parameter multipliers for the adjustments
    adj = new_setup['FiberSim_setup']['model']['manipulations']['adjustments']
    
    # Cycle through the adjustments setting the multiplier based on the
    # p_vector. We will make a new array of adjustments here to allow for
    # base variants
    
    new_adjustments = []
    p_counter = 0
    
    for (i, a) in enumerate(adj):
        
        span = a['factor_bounds'][1] - a['factor_bounds'][0]
        m = a['factor_bounds'][0] + p_vector[p_counter] * span
        a['multipliers'] = []
        if ('factor_mode' in a):
            if (a['factor_mode'] == "log"):
                a['multipliers'].append(np.power(10, m))
                del a['factor_mode']
        else:
            a['multipliers'].append(m)
        del a['factor_bounds']
        
        new_adjustments.append(a)
        
        p_counter = p_counter + 1
        
    # Now we need to check for base_variants
    # We will handle these by adding new elements to the multipliers list for
    # each adjustment
    if ('base_variants' in new_setup['FiberSim_setup']['model']['manipulations']):
        base_variants = new_setup['FiberSim_setup']['model']['manipulations']['base_variants']
        
        # Cycle through the base variants
        for (i,bv) in enumerate(base_variants):

            # Work out the index for the parameter we are adding
            new_mult_index = len(new_adjustments[0]['multipliers'])
            
            # And now the adjustments
            for (j, a) in enumerate(bv['adjustments']):

                # Work out what the new value will be
                if ('multipliers' in a):
                    new_value = a['multipliers'][0]
                else:
                    span = a['factor_bounds'][1] - a['factor_bounds'][0]
                    new_value = a['factor_bounds'][0] + p_vector[p_counter] * span
                    if (a['factor_mode'] == 'log'):
                        new_value = np.power(10, new_value)
                    
                    p_counter = p_counter + 1
                
                # Work out whether the adjustment is new or
                # matches an existing entry         
                matching_ind = return_matching_adjustment_index(
                                    new_adjustments, a)
                    
                # Branch depending on match
                if (matching_ind == -1):
                    # There's no match
                    # Duplicate the last multiplier for other adjustments
                    # Add in a 1, x muliplier for the test
                    
                    for na in new_adjustments:
                        y = na['multipliers']
                        if (len(y) <= new_mult_index):
                            na['multipliers'].append(y[-1])
                        else:
                            na['multipliers'][new_mult_index] = y[-1]
                        
                    # Now add in the new one
                    a['multipliers'] = [new_value]
                    a['multipliers'].insert(-1, 1)
                    
                    # Clean up the adjustment
                    if ('factor_bounds' in a):
                        del a['factor_bounds']
                        
                    if ('factor_mode' in a):
                        del a['factor_mode']
                    
                    new_adjustments.append(a)
                    
                else:
                    # There is a match.
                    # Add in fixed multiplier for the match
                    # Duplicate the last multiplier for the
                    # non-matching adjustments
                    
                    for (k,na) in enumerate(new_adjustments):
                        if (k == matching_ind):
                            # Match
                            y = na['multipliers']
                            if (len(y) <= new_mult_index):
                                na['multipliers'].append(new_value)
                            else:
                                na['multipliers'][new_mult_index] = new_value
                        else:
                            # Non-match
                            y = na['multipliers']
                            if (len(y) <= new_mult_index):
                                na['multipliers'].append(y[-1])
                            else:
                                na['multipliers'][new_mult_index] = \
                                    y[-1]
                                    
        # Delete the base variants that are no longer needed
        del new_setup['FiberSim_setup']['model']['manipulations']['base_variants']
            
    # Overwrite
    new_setup['FiberSim_setup']['model']['manipulations']['adjustments'] = \
       new_adjustments
        
    # Add Python objective call to characterization, with an
    # extra field for the progress file
    new_setup['FiberSim_setup']['characterization'][0] \
                        ['post_sim_Python_call'] = \
                    orig_setup['FiberSim_setup']['model']['fitting'] \
                       ['Python_objective_call']
                
    # Write the new setup to file
    file_name = json_analysis_file_string.split('/')[-1]
    
    new_setup_file_string = os.path.join(working_dir,
                                         file_name)
    new_setup_file_string = str(Path(new_setup_file_string).resolve().absolute())
    
    with open(new_setup_file_string, 'w') as f:
        json.dump(new_setup, f, indent=4)

    # Now run the characterize
    
    # First off, need to work out where FiberPy is
    FiberCpp_struct = new_setup['FiberSim_setup']['FiberCpp_exe']
    if (FiberCpp_struct['relative_to'] == 'this_file'):
        temp_dir = Path(json_analysis_file_string).parent.absolute()
    else:
        temp_dir = FiberCpp_struct['relative_to']
    FiberCpp_exe_file = os.path.join(temp_dir,
                                FiberCpp_struct['exe_file'])
    FiberCpp_exe_dir = Path(FiberCpp_exe_file)
    FiberPy_file = os.path.join(FiberCpp_exe_dir,
                                '../../code/FiberPy/FiberPy',
                                'FiberPy.py')
    FiberPy_file = str(Path(FiberPy_file).resolve().absolute())
    
    # Now characeterize
    characterize_model.characterize_model(new_setup_file_string)
    
    # test_function(p_vector)

    # At this point, error components should be in ../worker/trial_errors.xlsx
    trial_errors_file = os.path.join(working_dir, 'trial_errors.xlsx')
    trial_errors = pd.read_excel(trial_errors_file)
    
    # Make a dictionary from the trial_errors file
    prog_d = dict()
    prog_d['trial'] = progress_data['iteration']
    for i in range(len(p_vector)):
        prog_d['p_%i' % (i+1)] = p_vector[i]
    prog_d['error_total'] = trial_errors['error_total'][0]
    
    error_cpt_labels = trial_errors.columns
    for err_lab in error_cpt_labels:
        if ('error_cpt' in err_lab):
            prog_d[err_lab] = trial_errors[err_lab][0]
    
    # If there is a progress file, append the new entry to the dataframe
    trial_df = pd.DataFrame(data=prog_d, index=[0])
    if (os.path.isfile(progress_data['progress_file_string'])):
        prog_df = pd.read_excel(progress_data['progress_file_string'])
        prog_df = pd.concat([prog_df, trial_df])        
    else:
        # Make a new dataframe
        prog_df = trial_df
    
    # Clean the file and then write
    if (os.path.exists(progress_data['progress_file_string'])):
        os.remove(progress_data['progress_file_string'])
    prog_df.to_excel(progress_data['progress_file_string'], index=False)
            
    # Update
    progress_data['iteration'] = progress_data['iteration'] + 1
    
    # Take action if this is the best fit
    if (prog_d['error_total'] < progress_data['lowest_error']):
        progress_data['lowest_error'] = prog_d['error_total']
        update_best(progress_data,
                    trial_df,
                    json_analysis_file_string)

    # Copy files to current
    update_current(progress_data,
                   json_analysis_file_string)
        
    # Update a figure
    plot_progress(progress_data)
    
    # Return
    return prog_d['error_total']

def return_matching_adjustment_index(existing, test):
    """ Compares the test adjustment to the existing ones and
        returns an index for a match, or -1 if there is no match """
    
    for (i,a) in enumerate(existing):
        if (test['variable'] == a['variable']):
            if ('isotype' in a):
                # It's a kinetic parameter, check the other matches
                try:
                    if (test['isotype'] == a['isotype']) and \
                            (test['state'] == a['state']) and \
                            (test['transition'] == a['transition']) and \
                            (test['parameter_number'] == a['parameter_number']):
                        # It's a match
                        return i
                except:
                    pass
            else:
                # It's a non-kinetic parameter
                if (test['class'] == a['class']):
                    # It's a match
                    return i
    
    # No match
    return -1
    
def update_best(progress_data, trial_df, json_analysis_file_string):
    """ Updates best simulation """
    
    # Update the best dataframe
    if (os.path.isfile(progress_data['best_file_string'])):
        best_df = pd.read_excel(progress_data['best_file_string'])
        best_df = pd.concat([best_df, trial_df])
    else:
        best_df = trial_df
        
    # Clean the file
    if (os.path.exists(progress_data['best_file_string'])):
        os.remove(progress_data['best_file_string'])
        
    # Write it
    best_df.to_excel(progress_data['best_file_string'], index=False)
    
    # Now get the model file from the generated dir and copy that
    # To do this, we need to go back to the json analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        
    model_struct = json_data['FiberSim_setup']['model']
    if (model_struct['relative_to'] == 'this_file'):
        base_dir = str(Path(json_analysis_file_string).parent)
    else:
        base_dir = model_struct['relative_to']
    old_file = os.path.join(base_dir,
                            model_struct['fitting']['generated_folder'],
                            'model_1.json')
    new_file = os.path.join(progress_data['progress_folder'],
                            'best_model.json')
    shutil.copyfile(old_file, new_file)

    # Finally, get the files from the sim_output and copy them across
    char_struct = json_data['FiberSim_setup']['characterization'][0]
    if (char_struct['relative_to'] == 'this_file'):
        base_dir = str(Path(json_analysis_file_string).parent)
    else:
        base_dir = char_struct['relative_to']
    sim_output_dir = os.path.join(base_dir,
                                  char_struct['sim_folder'],
                                  'sim_output')
    sim_output_dir = str(Path(sim_output_dir).resolve().absolute())
    
    # Copy files across
    shutil.copytree(sim_output_dir, progress_data['progress_folder'],
                    dirs_exist_ok=True)
    
    # Call Python code if required
    fitting_struct = json_data['FiberSim_setup']['model']['fitting']
    if ('Python_best_call' in fitting_struct):
        model_struct = json_data['FiberSim_setup']['model']
        if (model_struct['relative_to'] == 'this_file'):
            base_dir = str(Path(json_analysis_file_string).parent)
        else:
            base_dir = model_struct['relative_to']
            
            
        python_cmd = 'python %s' % os.path.join(base_dir,
                                  fitting_struct['Python_best_call'])
        
        print(python_cmd)
        subprocess.call(python_cmd)
        
        
    
    # # Find the files and copy them to progress
    # for file in os.listdir(sim_output_dir):
    #     old_file = os.path.join(sim_output_dir, file)
    #     if (os.path.isfile(old_file)):
    #         new_file = os.path.join(progress_data['progress_folder'],
    #                                 file)
    #         shutil.copyfile(old_file, new_file)

def update_current(progress_data, json_analysis_file_string):
    """ Updates sim_output files for visualization """

    # Laod the json analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)

    # Get the files from the sim_output and copy them across
    char_struct = json_data['FiberSim_setup']['characterization'][0]
    if (char_struct['relative_to'] == 'this_file'):
        base_dir = str(Path(json_analysis_file_string).parent)
    else:
        base_dir = char_struct['relative_to']
    sim_output_dir = os.path.join(base_dir,
                                  char_struct['sim_folder'],
                                  'sim_output')
    sim_output_dir = str(Path(sim_output_dir).resolve().absolute())
    
    # Find the files and copy them to current
    for file in os.listdir(sim_output_dir):
        old_file = os.path.join(sim_output_dir, file)
        if (os.path.isfile(old_file)):
            new_file = os.path.join(progress_data['current_folder'],
                                    file)
            shutil.copyfile(old_file, new_file)
    
def plot_progress(progress_data):
    
    # Load in the progress data
    df = pd.read_excel(progress_data['progress_file_string'])
    
    # Make a 2 panel figure
    no_of_rows = 3
    no_of_cols = 1
    
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([5, 7])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
    
    # Plot progress
    ax[0].plot(df['trial'], np.log10(df['error_total']), 'b-')
    ax[0].set_xlabel('Iterations')
    ax[0].set_ylabel('log_{10} (Total Error)')
    
    # Plot component values
    best_row_ind = df['error_total'].idxmin()
    x = 1
    cur_label = []
    best_label = []
    for (i,col_lab) in enumerate(df.columns):
        if ('error_cpt' in col_lab):
            # Set the label
            if not (cur_label):
                cur_label = 'Current'
                best_label = 'Best'
            else:
                cur_label = '_Current'
                best_label = '_best'
            ax[1].plot(x, np.log10(df[col_lab].iloc[-1]), 'bs',
                       markerfacecolor='none',
                       label=cur_label)
            ax[1].plot(x, np.log10(df[col_lab].iloc[best_row_ind]), 'ro',
                       markerfacecolor='none',
                       label=best_label)
            x = x + 1
    ax[1].set_xlabel('Components')
    ax[1].set_ylabel('log_{10} (Component errors)')
    x_ticks = np.arange(x+1)
    ax[1].set_xticks(x_ticks)
    ax[1].legend(loc='upper left',
                 bbox_to_anchor=(1.0, 1.05),
                 fontsize='xx-small',
                 handlelength=1)
    
    # Plot parameter values
    x = 1
    cur_label = []
    best_label = []
    for (i,col_lab) in enumerate(df.columns):
        if ('p_' in col_lab):
            # Set the label
            if not (cur_label):
                cur_label = 'Current'
                best_label = 'Best'
            else:
                cur_label = '_Current'
                best_label = '_best'
                
            ax[2].plot(x, df[col_lab].iloc[-1], 'bs',
                       markerfacecolor = 'none',
                       label = cur_label)
            ax[2].plot(x, df[col_lab].iloc[best_row_ind], 'ro',
                       markerfacecolor = 'none',
                       label = best_label)
            x = x + 1
    ax[2].set_xlabel('Parameters')
    ax[2].set_ylabel('p value')
    y_ticks = [0, 1]
    ax[2].set_yticks(y_ticks)
    ax[2].legend(loc='upper left',
                 bbox_to_anchor=(1.0, 1.05),
                 fontsize='xx-small',
                 handlelength=1)
           
    
    # Save fig
    fig.savefig(os.path.join(progress_data['progress_folder'],
                             'progress.png'),
                bbox_inches='tight',
                dpi=200)

def test_function(p_vector):

    trial_errors_file = 'd:/ken/github/campbellmusclelab/models/fibersim/demo_files/fitting/ante/working/trial_errors.xlsx'
    
    g = np.asarray([0.1, 0.2, 0.3, 0.4, 0.5])
    
    
    d = dict()
    
    ev = np.NaN * np.ones(len(p_vector))
    
    rng = np.random.default_rng()
    
    for i in range(len(p_vector)):
        ev[i]= np.power(p_vector[i] - g[i], 2.0) + \
            (0.0 * (1+np.sin(100*p_vector[i]/3.14159))) + \
                (0.0*rng.random())
            
        d['error_cpt_%i' % (i+1)] = ev[i]
    d['error_total'] = np.sum(ev)
    
    print(d)
    
    df = pd.DataFrame(data=d, index=[0])
    
    df.to_excel(trial_errors_file, index=False)
    
            
    
    
    
    
    
    
    
    
    
    
    
        
    
    
        
    
        