# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 21:56:54 2024

@author: ken
"""

import os
import sys

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Add FiberSim code to path
sys.path.append('../../../code/fiberpy/fiberpy/package/modules/analysis')
import curve_fitting as cf

def return_fit(thread_id_dir):
    """ returns a single value defining the least squares fit
        between the current simulation and the target data """
       
    # Variables
    sim_data_folder = 'sim_data/sim_output'
    target_data_folder = '../../target_data/sim_output'
    
    pCa_data_file = 'pCa_analysis.xlsx'
    k_tr_data_file = 'k_tr_analysis.xlsx'
    
    trial_errors_file = 'working/trial_errors.xlsx'
    output_image_folder = 'sim_data/sim_output'
    
    pCa_eval_values = np.asarray([7, 6.3, 6.1, 5.9, 5.7, 5.5, 5.3, 4.5])
    
    k_tr_pCa_value = 4.5
    
    n_points_per_curve = 100
    
    # Code    

    # Adapt because files are relative to the thread_id_path
    sim_data_folder = Path(os.path.join(thread_id_dir, sim_data_folder)).resolve()
    sim_data_pCa_file = os.path.join(sim_data_folder, pCa_data_file)
    sim_data_k_tr_file = os.path.join(sim_data_folder, k_tr_data_file)
    
    target_data_folder = os.path.join(thread_id_dir, target_data_folder)
    target_data_pCa_file = os.path.join(target_data_folder, pCa_data_file)
    target_data_k_tr_file = os.path.join(target_data_folder, k_tr_data_file)
    
    trial_errors_file = os.path.join(thread_id_dir, trial_errors_file)

    # Open the sim data
    sim_data_pCa = pd.read_excel(sim_data_pCa_file)
    sim_data_k_tr = pd.read_excel(sim_data_k_tr_file)
    
    # Find the number of curves
    no_of_curves = sim_data_pCa['curve'].max()
    
    # Create a matrices with the interpolated pCa curves
    # Do this for both sim and target data
    x_fit = np.NaN * np.ones(n_points_per_curve)
    sim_curve = np.NaN * np.ones([n_points_per_curve, no_of_curves])
    target_curve = np.NaN * np.ones([n_points_per_curve, no_of_curves])
    
    n_eval_points = len(pCa_eval_values)
    sim_points = np.NaN * np.ones([n_eval_points, no_of_curves])
    target_points = np.NaN * np.ones([n_eval_points, no_of_curves])
    
    sim_k_tr = np.NaN * np.ones(no_of_curves)
    target_k_tr = np.NaN * np.ones(no_of_curves)
    
    # Fill them
    for ci in range(no_of_curves):
        
        # Get the sim data
        sim = sim_data_pCa[sim_data_pCa['curve'] == (ci+1)]
        pCa_raw = sim['hs_pCa'].to_numpy()
        f_raw = sim['hs_force'].to_numpy()
        
        # Fit
        fit_data = cf.fit_pCa_data(pCa_raw, f_raw, n_points_per_curve)
        
        # Fill in curve points
        if (ci == 0):
            x_fit = fit_data['x_fit']
            
        sim_curve[:, ci] = fit_data['y_fit']
        
        # Fill in the sim_pCa_values
        y_eval = y_pCa(pCa_eval_values,
                       fit_data['pCa_50'],
                       fit_data['n_H'],
                       fit_data['y_min'],
                       fit_data['y_amp'])
        
        sim_points[:, ci] = y_eval
        
        # Get the k_tr
        d_k_tr = sim_data_k_tr[(sim_data_k_tr['curve'] == (ci+1)) &
                               (sim_data_k_tr['pCa'] == k_tr_pCa_value)]
        
        sim_k_tr[ci] = d_k_tr['k_tr'].iloc[0]
        
        # Get the target data
        target_params = pd.read_excel(target_data_pCa_file,
                                      sheet_name = 'curve_%i' % (ci+1))

        # Fill the target curve and points
        target_curve[:, ci] = y_pCa(x_fit,
                                    target_params['pCa_50'].iloc[0],
                                    target_params['n_H'].iloc[0],
                                    target_params['y_min'].iloc[0],
                                    target_params['y_amp'].iloc[0])
        
        target_points[:, ci] = y_pCa(pCa_eval_values,
                                    target_params['pCa_50'].iloc[0],
                                    target_params['n_H'].iloc[0],
                                    target_params['y_min'].iloc[0],
                                    target_params['y_amp'].iloc[0])
        
        # Get the k_tr
        target_data_k_tr = pd.read_excel(target_data_k_tr_file)
        
        d_k_tr = target_data_k_tr[(target_data_k_tr['curve'] == (ci+1)) &
                               (target_data_k_tr['pCa'] == k_tr_pCa_value)]
        
        target_k_tr[ci] = d_k_tr['k_tr'].iloc[0]
                
    # Create a dict to store error results
    ed = dict()
    ed_total = 0
    
    # Calculate the least_squares_fit
    for ci in range(no_of_curves):
        
        y_target = target_points[:, ci]
        target_max = np.amax(y_target)
        
        y_sim = sim_points[:, ci]
        
        e = np.sum(np.power((y_target - y_sim) / target_max, 2.0))

        ed['error_cpt_%i' % (ci+1)] = e
        ed_total = ed_total + e
        
    # Add in the k_tr
    # Add in the least squares fit for the k_trs
    for ci in range(no_of_curves):
        e = np.power((sim_k_tr[ci] - target_k_tr[ci])/target_k_tr[ci], 2.0)
        ed['error_cpt_%i' % (no_of_curves + ci + 1)] = e
        ed_total = ed_total + e
        
    ed['error_total'] = ed_total

    # Make a dataframe
    df = pd.DataFrame(data=ed, index=[0])
    
    # Check the dir exists
    worker_parent_dir = Path(trial_errors_file).parent
    if not os.path.isdir(worker_parent_dir):
        os.makedirs(worker_parent_dir)
    
    # Clean the file and then write
    if (os.path.exists(trial_errors_file)):
        os.remove(trial_errors_file)
    df.to_excel(trial_errors_file, index=False)
    
    # Make a figure
    fig = plt.figure(constrained_layout=False)
    spec = fig.add_gridspec(nrows=1, ncols=1,
                      left=0.3, right=0.95, wspace=0.1,
                      bottom = 0.2)
    fig.set_size_inches([3.5, 2.5])
    ax = fig.add_subplot(spec[0, 0])
    
    for ci in range(no_of_curves):
        if (ci == 0):
            target_label = 'Expt'
            sim_label = 'Sim'
        else:
            target_label = '_'
            sim_label = '_'
            
        ax.plot(x_fit, target_curve[:, ci], 'k-', label=target_label)
        ax.plot(pCa_eval_values, target_points[:, ci], 'ko', label='_', mfc='none')
        ax.plot(x_fit, sim_curve[:, ci], 'b-', label=sim_label)
        ax.plot(pCa_eval_values, sim_points[:, ci], 'b^', label='_', mfc='none')
        
    ax.set_xlim(7, 4.5)
    ax.legend(fontsize=8)
    
    # Add k_tr
    y_lims = ax.get_ylim()
    y_anchor = y_lims[0] + (0.5 * (y_lims[1] - y_lims[0]))
    y_spacing = 0.15 * (y_lims[1] - y_lims[0])
    ax.text(6.85, y_anchor,
            'Expt k$_{tr}$ %.2f  %.2f' %
                (target_k_tr[0], target_k_tr[1]),
            fontsize=8,
            color = 'k')
    y_anchor = y_anchor - y_spacing
    ax.text(6.85, y_anchor,
            'Sim k$_{tr}$ %.2f %.2f' %
                (sim_k_tr[0], sim_k_tr[1]),
            fontsize=8,
            color = 'b')
    
    ofs = os.path.join(thread_id_dir, output_image_folder,
                        'curve_comparison.png')
    fig.savefig(ofs, dpi=200, bbox_inches='tight')


def y_pCa(x_data, pCa_50, n_H, y_min, y_amp):
    """ Returns a Hill curve from defined parameters """

    y = np.zeros(len(x_data))
    
    for (i, x) in enumerate(x_data):
        x = float(x)
        y[i] = y_min + \
            (y_amp * (np.power(np.power(10, -x), n_H)) /
             (np.power(np.power(10, -x), n_H) +
                      np.power(np.power(10, -pCa_50), n_H)))
            
    return y
    

if __name__ == "__main__":
    # Thread path gets included with function call as second argument    
    # Code needs this to work out where files are
    return_fit(sys.argv[1])
    