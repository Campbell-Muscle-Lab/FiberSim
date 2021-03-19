# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 16:51:53 2021

@author: ken
"""

import os
import json
import numpy as np
import pandas as pd
import math

from pathlib import Path

from shutil import copy

from collections import defaultdict

from scipy.optimize import minimize

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from package.modules.batch import batch as bat
from package.modules.analysis import curve_fitting as cv

class fitting():
    """ Class for fitting FiberSim to data """
    
    def __init__(self, optimization_json_file_string):
        """ Constructor for an optimization object """

        # Store the optimization_file_string
        self.optimization_json_file_string = optimization_json_file_string

        # Load the optimization task from file
        with open(self.optimization_json_file_string,'r') as f:
            json_data = json.load(f)

        # Create a dict with data for the optimization
        self.opt_data = dict()

        # Load optimization data
        if ('FiberSim_optimization' not in json_data):
            print('Error: FiberSim_optimization not found in %s' %
                  optimization_json_file_string)
            exit(1)

        opt_struct = json_data['FiberSim_optimization']
        for key in opt_struct.keys():
            if (not key == "files"):
                self.opt_data[key] = opt_struct[key]
            else:
                self.opt_data['files'] = dict()
                file_keys = opt_struct['files'].keys()
                for fk in file_keys:
                    if (not fk == 'relative_to'):
                        fs = opt_struct['files'][fk]
                        if (not opt_struct['files']['relative_to']):
                            fs = os.path.abspath(fs)
                        elif (opt_struct['files']['relative_to'] ==
                              'this_file'):
                            base_directory = \
                                Path(optimization_json_file_string).\
                                    parent.absolute()
                            fs = os.path.join(base_directory, fs)
                        else:
                            base_directory = opt_struct['files']['relative_to']
                            fs = os.path.join(base_directory, fs)
                        self.opt_data['files'][fk] = fs

        # Load batch structure
        if ('FiberSim_batch' not in json_data):
            print('Error: FiberSim_batch not founs in %s' %
                  optimization_json_file_string)
            exit(1)
                        
        self.batch_structure = json_data['FiberSim_batch']
        # Update batch structure for internal use
        job_data = self.batch_structure['job']
        new_job_data = []
        for i, j in enumerate(job_data):
            d = dict()
            for f in ['model_file', 'options_file',
                  'protocol_file', 'results_file']:
                fs = j[f]
                if (not j['relative_to']):
                    fs = os.path.abspath(fs)
                elif (j['relative_to'] == 'this_file'):
                    base_directory = Path(optimization_json_file_string). \
                        parent.absolute()
                    fs = os.path.join(base_directory, fs)
                else:
                    base_directory = j['relative_to']
                    fs = os.path.join(base_directory, fs)
                d[f] = fs
            new_job_data.append(d)
        self.batch_structure['job'] = new_job_data

        # Load length initial conditions if specified
        if('initial_delta_hsl' in json_data):
            self.initial_delta_hsl = json_data['initial_delta_hsl'] 

        # Load constraint structure if specified
        if('constraint' in json_data):
            self.constraint = json_data['constraint'] 

        # Load parameter data
        if ('parameter' not in json_data):
            print('Error: no parameter specified in %s' %
                  optimization_json_file_string)
            exit(1)

        # Create data for optimization
        par_data = json_data['parameter']
        self.p_data = []
        self.p_vector = []
        for p in par_data:
            p_obj = parameter(p)
            self.p_data.append(p_obj)
            self.p_vector.append(p_obj.data['p_value'])

        # Add in contraints 
        if ('constraint' in json_data):
            for constr in json_data['constraint']:
                if('parameter_multiplier' in constr):
                    for param_mult in constr["parameter_multiplier"]:
                        p_obj = parameter(param_mult)
                        p_obj.data['name'] = "multi_" + p_obj.data['name']
                        self.p_data.append(p_obj)
                        self.p_vector.append(p_obj.data['p_value'])

        self.p_vector = np.array(self.p_vector)
        self.global_fit_values = np.array(0)
        self.best_fit_value = []
        self.best_fit_data = []
        self.best_p_vector = []
        
        # Start controller
        self.fit_controller()

    def fit_controller(self):
        """ Controls fitting routines """

        print('Initialising fit controller')
        self.global_fit_values = []
        self.best_fit_value = np.inf

        if (self.opt_data['optimizer'] == 'particle_swarm'):
            # Set up pyswarms
            import pyswarms as ps
            options = {'c1': 0.5, 'c2': 0.5, 'w': 0.5}
            max_bounds = np.ones(len(self.p_vector))
            min_bounds = np.zeros(len(self.p_vector))
            bounds = (min_bounds, max_bounds)
            print(bounds)
            optimizer = ps.single.GlobalBestPSO(
                n_particles = round(1.5 * len(self.p_vector)),
                dimensions = len(self.p_vector),
                options=options,
                bounds = bounds)
            # Initialise
            cost, pos = optimizer.optimize(self.pso_wrapper,
                                           iters = np.round(100*len(self.p_vector)))
        else:
            res = minimize(self.fit_worker, self.p_vector,
                            method='Nelder-Mead',
                            tol=1e-3)

    def pso_wrapper(self, p_array):
        """ Runs an iteration from PySwarms """

        # Code
        fit_error = []
        for row in p_array:
            e = self.fit_worker(row)
            fit_error.append(e)

        # Return
        return fit_error


    def fit_worker(self, p_vector):
        """ Runs a batch process, evaluates fits
        and updates progress """
        
        print(p_vector)

        # First update the p_data
        for i, p in enumerate(self.p_data):
            p.data['p_value'] = p_vector[i]
            p.data['value'] = p.return_parameter_value()

        # Now update the model worker files
        for i, j in enumerate(self.batch_structure['job']):
            self.update_model_worker_file(i, j)

        # Now run batch
        print('Running %.0f jobs as a batch process' %
              np.size(self.batch_structure['job']))
        bat.run_batch(self.optimization_json_file_string)

        # Finally evaluate the fits and append
        print('Evaluating fits for the batch process')
        fit_data = self.evaluate_fit()
        fit_error = fit_data['fit_error']
        self.global_fit_values.append(fit_error)

        # Check for best error and implement
        if (fit_error <= self.best_fit_value):
            self.best_fit_value = fit_error
            self.best_fit_data = fit_data
            self.best_p_vector = p_vector

            # Update the best model files
            if not os.path.exists(self.opt_data['files']['best_model_folder']):
                os.makedirs(self.opt_data['files']['best_model_folder'])
            for i, j in enumerate(self.batch_structure['job']):
                ofs = os.path.split(j['model_file'])[-1]
                nfs = os.path.join(
                    self.opt_data['files']['best_model_folder'], str(i+1))
                if not os.path.exists(nfs):
                    os.makedirs(nfs)
                print('Copying model file from\n%s\nto\n%s' %
                      (j['model_file'], nfs))
                copy(j['model_file'], nfs)

            # Save the best opt file
            json_data = dict()
            json_data['FiberSim_optimization']  = self.opt_data;
            json_data['FiberSim_batch'] = self.batch_structure
            json_data['parameter']=[]
            for p in self.p_data:
                json_data['parameter'].append(p.data)

            dir_name = os.path.dirname(os.path.abspath(
                self.opt_data['files']['best_opt_file']))
            if (not os.path.isdir(dir_name)):
                os.makedirs(dir_name)
            with open(self.opt_data['files']['best_opt_file'], 'w') as f:
                json.dump(json_data, f, indent=4)

        # Save figures
        if ('figure_current_fit_file' in self.opt_data['files']):
            self.create_figure_current_fit(fit_data)
        if ('figure_fit_progress_file' in self.opt_data['files']):
            self.create_figure_fit_progress()

        # If single run, abort
        if (self.opt_data['single_run'] == 1):
            print('Single run set in optimization file - now exiting')
            exit(1)

        # Return fit
        return fit_error


    def evaluate_fit(self):
        """ Branches to different fit functions depending on mode """
        
        if (self.opt_data['fit_mode'] == 'fit_in_time_domain'):
            fit_data = self.evaluate_time_domain_fit()

        if (self.opt_data['fit_mode'] == 'fit_pCa_curve'):
            fit_data = self.evaluate_pCa_curve_fit()

        return fit_data


    def evaluate_time_domain_fit(self):
        """ Evaluates time domain fit """

        # Deduce the number of jobs
        no_of_jobs = np.size(self.batch_structure['job'])
        
        # Create a dictionary to hold data describing the fit
        fit_data = dict()
        fit_data['job'] = []
        fit_data['job_errors'] = np.zeros(no_of_jobs)
        fit_data['fit_error'] = []

        # Loop through jobs
        for i, j in enumerate(self.batch_structure['job']):
            sim_file_string = os.path.join(j['output_folder'], 'results.txt')
            d = pd.read_csv(sim_file_string, sep='\t')
            n_data = np.size(d[self.opt_data['fit_variable']])
            # Load in target and add to end of data structure
            target_file_string = j['target_file_string']
            target = pd.read_csv(target_file_string)
            d['target'] = np.NaN
            # Now copy target into data structure
            n_target = np.size(target)
            vi = np.arange(n_data-n_target, n_data, 1).astype('int')
            y_target = target.to_numpy(copy=True)
            max_target = np.amax(np.abs(y_target))
            d.loc[vi, 'target'] = y_target
            # Now calculate the sum of square differnce
            y_dif = (d.loc[vi, 'target'] - d.loc[vi, self.opt_data['fit_variable']]) / \
                        max_target
            fit_data['job_errors'][i] = \
                np.sqrt(np.sum(np.power(y_dif, 2))) / n_target

            # Keep time, force and target data and add to
            # fit_data
            fit_data['job'].append(
                d.filter(items=['time', 'target', self.opt_data['fit_variable']]))

        # Global error
        fit_data['fit_error'] = np.sum(fit_data['job_errors'])

        # Return results
        return fit_data

    def evaluate_pCa_curve_fit(self):
        """ Evaluates pCa curve fit """

        # Get the target data

        target_file_string = self.opt_data['files']['target_file']
        target = pd.read_excel(target_file_string, engine='openpyxl')
        
        # Deduce the number of jobs
        # no_of_jobs = np.size(self.batch_structure['job'])
        
        # Create a dictionary to hold data describing the fit
        fit_data = dict()
        fit_data['job'] = []
        fit_data['fit_error'] = []

        # Create dictionnaries for pCa, target data and calculated fit_variable

        target_data = defaultdict(list)
        calculated_data = defaultdict(list)
        pCa_data = defaultdict(list)

        # Create dictionnaries for the errors calculation

        max_y_target = defaultdict(list)
        y_dif = defaultdict(list)

        # Loop through the pCa curves
        for i, j in enumerate(target['curve']):

            # Check that each cell contains a numerical value
            if math.isnan(j): 
                break
            
            # Get pCa values
            pCa_data[j-1].append(target['pCa'][i])
            
            # Get target data
            target_data[j-1].append(target[self.opt_data['fit_variable']][i])

            # Get max element of target data
            max_y_target[j-1] = np.amax(np.abs(target_data[j-1]))
                        
            # Get simulation data
            job_name = self.batch_structure['job'][i]
            sim_file_string = job_name['results_file']
            d = pd.read_csv(sim_file_string, sep='\t')

            # Get last element of fit_variable
            calculated_data[j-1].append(d[self.opt_data['fit_variable']].iloc[-1])

        no_of_curves = len(max_y_target)
        n_data = len(calculated_data[0])
        
        fit_data['job_errors'] = np.zeros(no_of_curves)

        # Calculate error for each pCa curve
        for i, max_y in enumerate(max_y_target):       
            
            y_dif[i] = (np.array(calculated_data[i]) - 
                        np.array(target_data[i]))/np.array(max_y_target[i]) 

            fit_data['job_errors'][i] = \
                np.sqrt(np.sum(np.power(y_dif[i], 2))) / np.size(y_dif[i])

            df = pd.DataFrame()
            df['calculated_data'] = calculated_data[i]
            df['target_data'] = target_data[i]
            df['pCa'] = pCa_data[i]
            # Keep curve, pCa, target and calculated data and add to
            # fit_data
            fit_data['job'].append(df)

        # Global error
        fit_data['fit_error'] = np.sum(fit_data['job_errors'])

        # Return results
        return fit_data


    def create_figure_fit_progress(self):
        """ Creates a figure showing the progress of the fit """

        fig = plt.figure(constrained_layout=True)
        fig.set_size_inches([7, 7])
        spec = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
        ax = []
        ax.append(fig.add_subplot(spec[0, 0]))

        x = np.arange(0, np.size(self.global_fit_values), 1)
        if (np.size(x)==1):
            marker = 'bo'
        else:
            marker = 'b-'
        ax[0].plot(x, np.log(self.global_fit_values), marker)
        ax[0].set_xlabel('Fit values')
        ax[0].set_ylabel('log$_{10}$\nfit\nerror')

        # Save figure
        print('Saving fit progress to %s' %
               self.opt_data['files']['figure_fit_progress_file'])
        # Check folder exists and make it if not
        dir_name = os.path.dirname(os.path.abspath(
            self.opt_data['files']['figure_fit_progress_file']))
        if (not os.path.isdir(dir_name)):
                os.makedirs(dir_name)
        fig.savefig(self.opt_data['files']['figure_fit_progress_file'])
        plt.close()


    def create_figure_current_fit(self, fit_data):

        fig = plt.figure(constrained_layout=True)
        fig.set_size_inches([7, 7])
        spec = gridspec.GridSpec(nrows=3, ncols=1, figure=fig)
        spec.update(left=0.3, right=0.7)
        ax=[]
        if (self.opt_data['fit_mode'] == 'fit_in_time_domain'):
            # Fit trace against target plotted v time
            ax.append(fig.add_subplot(spec[0,0]))
            
            for j in fit_data['job']:
                ax[0].plot(j['time'], j['target'], 'k-')
                ax[0].plot(j['time'], j[self.opt_data['fit_variable']], 'b-')
            # Add in best_fit
            if (self.best_fit_data):
                for j in self.best_fit_data['job']:
                    ax[0].plot(j['time'], j[self.opt_data['fit_variable']], 'r-')    

        if (self.opt_data['fit_mode'] == 'fit_pCa_curve'):
            # Fit trace against target pCa curve
            ax.append(fig.add_subplot(spec[0,0]))

            cf=[]
            cf2=[]
            for j in fit_data['job']:
                ax[0].plot(j['pCa'], j['target_data'], 'ko')
                # Add in curve_fitting
                res = cv.fit_pCa_data(j['pCa'],j['target_data'])
                ax[0].plot(res["x_fit"], res["y_fit"], 'k-')
                cf.append(res)

                ax[0].plot(j['pCa'], j['calculated_data'], 'bo')
                # Add in curve_fitting
                res = cv.fit_pCa_data(j['pCa'],j['calculated_data'])
                ax[0].plot(res["x_fit"], res["y_fit"], 'b-')
                cf2.append(res)

            # Add in best_fit
            if (self.best_fit_data):
                cf3=[]
                for j in self.best_fit_data['job']:
                    ax[0].plot(j['pCa'], j['calculated_data'], 'ro')

                    # Add in curve_fitting
                    res = cv.fit_pCa_data(j['pCa'],j['calculated_data'])
                    # #x_data = np.power(10,-res["x_fit"])
                    ax[0].plot(res["x_fit"], res["y_fit"], 'r-')
                    cf3.append(res)

            ax[0].invert_xaxis()
            ax[0].set_xlim((7,4))
            ylim = ax[0].get_ylim()
            y_anchor = 0.9
            y_spacing = 0.1
            for c in cf:
                ax[0].text(6.9, y_anchor * ylim[1],
                       ('pCa50: %.2f n_H: %.2f' % (c['pCa_50'], c['n_H'])),
                       color='black')
                y_anchor = y_anchor - y_spacing
            for c in cf2:
                ax[0].text(6.9, y_anchor * ylim[1],
                       ('pCa50: %.2f n_H: %.2f' % (c['pCa_50'], c['n_H'])),
                       color='blue')
                y_anchor = y_anchor - y_spacing

            if (self.best_fit_data):
                for c in cf3:
                    ax[0].text(6.9, y_anchor * ylim[1],
                       ('pCa50: %.2f n_H: %.2f' % (c['pCa_50'], c['n_H'])),
                       color='red')
                    y_anchor = y_anchor - y_spacing

        ax.append(fig.add_subplot(spec[1,0]))
        for i, j in enumerate(fit_data['job_errors']):
            ax[1].plot(i+1, np.log10(j), 'ko')
        ax[1].set_xlabel('Job number')
        ax[1].set_ylabel('log$_{10}$\nJob\nerror')

        # Plot the p vector
        ax.append(fig.add_subplot(spec[2,0]))
        for i, p in enumerate(self.p_data):
            y = np.size(self.p_data) - i
            ax[2].plot(p.data['p_value'], y, 'bo')
            ax[2].plot(self.best_p_vector[i], y, 'rs')
            name_string = p.data['name']
            if (name_string == 'm_kinetics'):
                if ('extension' in p.data):
                    name_string = name_string + \
                        ('\nstate_%i_extension' % p.data['state'])
                if ('parameter_index' in p.data):
                    name_string = name_string + \
                        ('\ntransition_%i_to_%i[%i]' %
                         (p.data['old_state'], p.data['new_state'],
                          p.data['parameter_index']))
            ax[2].text(-2, y, name_string,
                       horizontalalignment = 'left',
                       clip_on=False)
            ax[2].text(3, y, ('%4g' % p.data['value']),
                       horizontalalignment = 'right',
                       clip_on=False)
        for i in np.arange(-1, 3, 1):
            ax[2].plot([i, i], [1, 1 + np.size(p.data)],
                       'k:')

        # Save figure
        print('Saving current fit to %s' %
               self.opt_data['files']['figure_current_fit_file'])
        # Check folder exists and make it if not
        dir_name = os.path.dirname(os.path.abspath(
            self.opt_data['files']['figure_current_fit_file']))
        if (not os.path.isdir(dir_name)):
                os.makedirs(dir_name)
        fig.savefig(self.opt_data['files']['figure_current_fit_file'])
        plt.close()


    def update_model_worker_file(self, job_numb, job_data):
        # Writes a new model worker file based on the p vector
        
        # First load in the model_template
        with open(self.opt_data['files']['model_template_file'], 'r') as f:
            model_template = json.load(f)

        # Nested loop through jobs and parameters
        for i, p in enumerate(self.p_data):
            new_value = p.return_parameter_value()
            if (p.data['name'] == 'm_kinetics'):
                self.replace_m_kinetics(model_template, p.data, new_value)
            elif (p.data['name'] == 'c_kinetics'):
                print('replace_c_kinetics')
            else:
                model_template = self.replace_item(model_template,
                                                   p.data['name'], new_value)

        # Check for length conditions 
        if(hasattr(self, 'initial_delta_hsl')):
            model_template["muscle"]['initial_hs_length'] = \
                model_template["muscle"]['initial_hs_length'] + self.initial_delta_hsl[job_numb]

        # Check contraints
        if(hasattr(self, 'constraint')):

            for constr in self.constraint:
                # Check if constraint is associated with the job model currently being updated
                if constr["job_number"] == job_numb + 1: # first job has job_numb = 0
                    # Check for parameter multiplier condition
                    if "parameter_multiplier" in constr:
                        for multi_data in constr["parameter_multiplier"]:
                            # find model file associated with base_job_number
                            base_job = self.batch_structure["job"][multi_data["base_job_number"]-1]
                            with open(base_job['model_file'], 'r') as f:
                                model_base = json.load(f)

                            # get base parameter value from model_base
                            base_model, base_value = self.find_item(model_base, multi_data['name'], 0) 
                            # get multiplier value
                            for i, p in enumerate(self.p_data):
                                if ("multi_" + multi_data['name']) in p.data['name']:
                                    multiplier_value = p.return_parameter_value()  
                            new_value = base_value * multiplier_value 
                            # replace new parameter value 
                            model_template = self.replace_item(model_template, multi_data['name'], new_value)
                    # Check for parameter copy condition
                    if "parameter_copy" in constr:
                        for copy_data in constr["parameter_copy"]:
                            # find model file associated with copy_job_number
                            copy_job = self.batch_structure["job"][copy_data["copy_job_number"]-1]

                            with open(base_job['model_file'], 'r') as f:
                                model_base = json.load(f)
                            # get base parameter value from model_base
                            copy_model, new_value = self.find_item(model_base, copy_data['name'], 0) # take the param value 

                            # replace new parameter value  
                            model_template = self.replace_item(model_template, copy_data['name'], new_value)
               
        # Now write updated model to file
        model_fs = job_data['model_file']

        # Check the folder exists and make it if required
        dir_name = os.path.dirname(model_fs)
        
        print('dir_name %s' % dir_name)
        if (not os.path.isdir(dir_name)):
                os.makedirs(dir_name)
        # Now dump file
        with open(model_fs,'w') as f:
            json.dump(model_template, f, indent=4)
            print('file_written')

    def replace_m_kinetics(self, model_template, par_struct, new_value):
        """ replace m_kinetics transition element """
        m_kinetics = model_template['m_kinetics']
        sch = m_kinetics['scheme']

        # Handle extension first
        if ('extension' in m_kinetics):
            s_new = []
            for i, s in enumerate(sch):
                if (s['number'] == par_struct['state']):
                    s['extension'] = new_value
                s_new.append(s)
            m_kinetics['scheme'] = s_new
        else:
        # Now do rate function parameters
            s_new = []
            for i, s in enumerate(sch):
                t_new = []
                for j, t in enumerate(s['transition']):
                    if ((s['number'] == par_struct['old_state']) and
                        (t['new_state'] == par_struct['new_state'])):
                            t['rate_parameters'][par_struct['parameter_index']] = \
                                new_value
                    t_new.append(t)
                s_new.append(s)
            m_kinetics['scheme'] = s_new

        # Return
        model_template['m_kinetics'] = m_kinetics
        return model_template

    def replace_item(self, obj, key, replace_value):
        # See https://stackoverflow.com/questions/45335445/recursively-replace-dictionary-values-with-matching-key
        for k, v in obj.items():
            if isinstance(v, dict):
                obj[k] = self.replace_item(v, key, replace_value)
        if key in obj:
            obj[key] = replace_value
        return obj

    def find_item(self, obj, key, val):
        # find value (val) associated with a given key (key) in a nested structure (obj)
        for k, v in obj.items():
            if isinstance(v, dict):
                obj[k], val = self.find_item(v, key, val)
        if key in obj:
            val = obj[key] 
        return obj, val


class parameter():
    """ Class for a parameter value """

    def __init__(self, parameter_struct):
        """ Constructor a parameter object """

        self.data = dict()
        for key in parameter_struct.keys():
            self.data[key] = parameter_struct[key]
        self.data['value'] = self.return_parameter_value()

    def return_parameter_value(self):
        """ Function returns the value of the parameter """

        temp_value = np.mod(self.data['p_value'], 2)
        if (temp_value < 1):
            parameter_value = self.data['min_value'] + \
                temp_value * (self.data['max_value'] - self.data['min_value'])
        else:
            parameter_value = self.data['max_value'] - \
                (temp_value - 1) * \
                    (self.data['max_value'] - self.data['min_value'])

        if (self.data['p_mode'] == 'log'):
            parameter_value = np.power(10, parameter_value)

        return parameter_value

