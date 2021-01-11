# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 16:51:53 2021

@author: ken
"""

import os
import json
import numpy as np
import pandas as pd

from shutil import copyfile

from scipy.optimize import minimize

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from package.modules.batch import batch as bat

class fitting():
    """ Class for fitting FiberSim to data """
    
    def __init__(self, optimization_json_file_string):
        """ Constructor for an optimization object """

        # Load the optimization task from file
        with open(optimization_json_file_string,'r') as f:
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
            self.opt_data[key] = opt_struct[key]
        
        # Load batch structure
        if ('FiberSim_batch' not in json_data):
            print('Error: FiberSim_batch not foudn in %s' %
                  optimization_json_file_string)
            exit(1)
                        
        self.batch_structure = json_data['FiberSim_batch']

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
        self.p_vector = np.array(self.p_vector)
        self.global_fit_values = np.array(0)
        self.best_fit_value = []
        self.best_fit_data = []
        self.best_p_vector = []



    def fit_controller(self):
        """ Controls fitting routines """

        print('Initialising fit controller')
        self.global_fit_values = []
        self.best_fit_value = np.inf
        
        res = minimize(self.fit_worker, self.p_vector,
                       method='Nelder-Mead',
                       tol=1e-3)


    def fit_worker(self, p_vector):
        """ Runs a batch process, evaluates fits
        and updates progress """

        # First update the p_data
        for i, p in enumerate(self.p_data):
            p.data['p_value'] = p_vector[i]

        # Now update the model worker files
        for j in self.batch_structure['job']:
            self.update_model_worker_file(j)

        # Now run batch
        print('Running %.0f jobs as a batch process' %
              np.size(self.batch_structure['job']))
        bat.run_batch(batch_structure = self.batch_structure)

        # Finally evaluate the fits and append
        print('Evaluating fits for the batch process')
        fit_data = self.evaluate_fit()
        fit_error = fit_data['fit_error']
        self.global_fit_values.append(fit_error)

        # Chegck for best error and implement
        if (fit_error <= self.best_fit_value):
            self.best_fit_value = fit_error
            self.best_fit_data = fit_data
            self.best_p_vector = p_vector

            # Update the best model files
            if not os.path.exists(self.opt_data['best_model_folder']):
                os.makedirs(self.opt_data['best_model_folder'])
            for j in self.batch_structure['job']:
                ofs = os.path.split(j['model_file_string'])[-1]
                nfs = os.path.join(
                    self.opt_data['best_model_folder'], ofs)
                print('Copying model file from\n%s\nto\n%s' %
                      (j['model_file_string'], nfs))
                copyfile(j['model_file_string'], nfs)

            # Save the best opt file
            json_data = dict()
            json_data['FiberSim_optimization']  = self.opt_data;
            json_data['FiberSim_batch'] = self.batch_structure
            json_data['parameter']=[]
            for p in self.p_data:
                json_data['parameter'].append(p.data)
            dir_name = os.path.dirname(os.path.abspath(
                self.opt_data['best_opt_file_string']))
            if (not os.path.isdir(dir_name)):
                os.makedirs(dir_name)
            with open(self.opt_data['best_opt_file_string'], 'w') as f:
                json.dump(json_data, f, indent=4)

        # Save figures
        if ('figure_current_fit' in self.opt_data):
            self.create_figure_current_fit(fit_data)
        if ('figure_fit_progress' in self.opt_data):
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

        return fit_data


    def evaluate_time_domain_fit(self):
        """ Evaluates time domain fit """

        # Deduce the number of jobs
        no_of_jobs = np.size(self.batch_structure['job'])
        
        # Create a dictionary to hold data decribing the fit
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
               self.opt_data['figure_fit_progress'])
        fig.savefig(self.opt_data['figure_fit_progress'])


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
            ax[2].text(-2, y, p.data['name'],
                       horizontalalignment = 'left',
                       clip_on=False)
            ax[2].text(3, y, ('%4g' % p.data['p_value']),
                       horizontalalignment = 'right',
                       clip_on=False)
        for i in np.arange(-1, 3, 1):
            ax[2].plot([i, i], [1, 1 + np.size(p.data)],
                       'k:')
        
        # Save figure
        print('Saving current fit to %s' %
               self.opt_data['figure_current_fit'])
        fig.savefig(self.opt_data['figure_current_fit'])


    def update_model_worker_file(self, job_data):
        # Writes a new model worker file based on the p vector
        
        # First load in the model_template
        with open(self.opt_data['model_template_file_string'], 'r') as f:
            model_template = json.load(f)

        # Nested loop through jobs and parameters
        for i, p in enumerate(self.p_data):
            new_value = p.return_parameter_value()
            model_template = self.replace_item(model_template,
                              p.data['name'], new_value)

        # Now write updated model to file
        with open(job_data['model_file_string'],'w') as f:
            json.dump(model_template, f, indent=4)


    def replace_item(self, obj, key, replace_value):
        # See https://stackoverflow.com/questions/45335445/recursively-replace-dictionary-values-with-matching-key
        for k, v in obj.items():
            if isinstance(v, dict):
                obj[k] = self.replace_item(v, key, replace_value)
        if key in obj:
            obj[key] = replace_value
        return obj


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

