# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:41:41 2023

@author: Campbell
"""

import os
os.environ["OMP_NUM_THREADS"] = '1'
import json
import sys
import copy
import re
import shutil

from sklearn.cluster import KMeans

import numpy as np
import pandas as pd

import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from natsort import natsorted

from pathlib import Path

def pair_plot(excel_data_file_string, analysis_folder):
    
    d = pd.read_excel(excel_data_file_string)
    
    sns.set_theme(style = 'ticks')
    sns_plot = sns.pairplot(d)
    
    fig = sns_plot.fig
    
    ofs = os.path.join(analysis_folder, 'pair_plot.png')
    print('Saving pair_plot_figure to: %s' % ofs)
    fig.savefig(ofs, bbox_inches='tight')

def analyze_simulations(top_data_folder, sim_file_string, analysis_folder):
    """ Runs through simulations, pulling off summary statistics
        for each one, returns results in a pandas dataframe """
       
    # Pull off the data directories
    data_dirs = [f for f in Path(top_data_folder).iterdir() if f.is_dir()]
    for (i,d) in enumerate(data_dirs):
            
        dfs = os.path.join(os.path.abspath(str(d)),
                       sim_file_string)
        
        sample_results = analyze_simulation(dfs, i)
        
        # Store the data
        if (i==0):
            collated_results = sample_results
        else:
            collated_results = pd.concat([collated_results, sample_results],
                                         ignore_index = True)
            
    # Write results to file
    results_file_string = os.path.join(analysis_folder, 'analysis.xlsx')
    collated_results.to_excel(results_file_string, index=False)
        
            
def analyze_simulation(data_file_string,
                       sim_index,
                       rolling_n = 20):
    """ Analyze a single simulation, returns results in a pandas dataframe """
    
    # Load the data
    d = pd.read_csv(data_file_string, sep='\t', na_values=['-nan(ind'])
    
    # Make a figure
    # Set up a figure
    no_of_rows = 5
    no_of_cols = 1
        
    Ca_row = 1
    dCadt_row = 2
    hsl_row = 3
    dhsldt_row = 4
    d2Cadt2_row = 5
    
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([3, 10])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
            
    # Make the additional columns
    d['Ca'] = np.power(10, -d['hs_1_pCa'])
    d['dCadt'] = np.gradient(d['Ca'], d['time'])
    d['dCadt'] = d['dCadt'].rolling(rolling_n).sum(). \
        fillna(method='bfill').fillna(method='ffill')
    d['d2Cadt2'] = np.gradient(d['dCadt'], d['time'])
    
    # Get some values
    res = dict()
    res['Ca_max'] = d['Ca'].max()
    res['Ca_max_time_s'] = d['time'].iloc[d['Ca'].argmax()]
    res['dCadt_max'] = d['dCadt'].max()
    res['dCadt_max_time_s'] = d['time'].iloc[d['dCadt'].argmax()]
    res['dCadt_min'] = d['dCadt'].min()
    res['dCadt_min_time_s'] = d['time'].iloc[d['dCadt'].argmin()]
    
    # Cluster the dCadt
    d2 = d[(d['time'] <= res['Ca_max_time_s'])]
    dCadt = d2['dCadt'].values
    dCadt = dCadt.reshape(-1,1)
    kmeans = KMeans(n_clusters = 2, n_init=5).fit(dCadt)
    
    # Get the value of the last point (probably in the peak)
    trans_label = kmeans.labels_[-1]
    
    vi = np.nonzero(kmeans.labels_ == trans_label)[0]
    vi = vi[0] - 1
   
    # vi = np.nonzero((d['d2Cadt2'] > (0.01 * d['d2Cadt2'].max())) &
    #                 (d['time'] < res['Ca_max_time_s']))[0]
    # vi = vi[0]    
    res['stim_time_s'] = d['time'].iloc[vi]
    res['Ca_initial'] = d['Ca'].iloc[vi]
    
    # Now look at hs_length
    res['hs_length_initial'] = d['hs_1_length'].iloc[vi]
    # We have to be cunning to get the min length after the stimulus
    d['temp'] = d['hs_1_length']
    d.loc[d['time'] < res['stim_time_s'], 'temp'] = np.NaN
    
    res['hs_length_min'] = np.nanmin(d['temp'].values)
    res['hs_length_min_time_s'] = d['time'].iloc[np.nanargmin(d['temp'].values)]
    
    # And the times to relative values
    
    thresh = [0.1, 0.5, 0.9]
    for t in thresh:
        lab = ('Ca_up_%.0f' % (100 * t))
        y_thresh = res['Ca_initial'] + (t * (res['Ca_max'] - res['Ca_initial']))
        vi = last_point_before(d['time'], d['Ca'],
                               res['stim_time_s'], res['Ca_max_time_s'],
                               y_thresh, 'below')
        if (vi >= 0):
            res['%s_time_s' % lab] = d['time'].iloc[vi]
            res['%s' % lab] = d['Ca'].iloc[vi]
        else:
            res['%s_time_s' % lab] = np.NaN
            res['%s' % lab] = np.NaN
 
    for t in thresh:
        lab = ('Ca_down_%.0f' % (100 * t))
        y_thresh = res['Ca_initial'] + (t * (res['Ca_max'] - res['Ca_initial']))
        vi = last_point_before(d['time'], d['Ca'],
                               res['Ca_max_time_s'], d['time'].max(),
                               y_thresh, 'above')
        if (vi >= 0):
            res['%s_time_s' % lab] = d['time'].iloc[vi]
            res['%s' % lab] = d['Ca'].iloc[vi]
        else:
            res['%s_time_s' % lab] = np.NaN
            res['%s' % lab] = np.NaN
            
    for t in thresh:
        lab = ('hs_length_down_%.0f' % (100 * t))
        y_thresh = res['hs_length_initial'] + \
                    (t * (res['hs_length_min'] - res['hs_length_initial']))
        vi = last_point_before(d['time'], d['hs_1_length'],
                                res['stim_time_s'], res['hs_length_min_time_s'],
                                y_thresh, 'above')
        if (vi >= 0):
            res['%s_time_s' % lab] = d['time'].iloc[vi]
            res['%s' % lab] = d['hs_1_length'].iloc[vi]
        else:
            res['%s_time_s' % lab] = np.NaN
            res['%s' % lab] = np.NaN
            
    for t in thresh:
        lab = ('hs_length_up_%.0f' % (100 * t))
        y_thresh = res['hs_length_initial'] + \
                    (t * (res['hs_length_min'] - res['hs_length_initial']))
        vi = last_point_before(d['time'], d['hs_1_length'],
                                res['hs_length_min_time_s'], d['time'].max(),
                                y_thresh, 'below')
        if (vi >= 0):
            res['%s_time_s' % lab] = d['time'].iloc[vi]
            res['%s' % lab] = d['hs_1_length'].iloc[vi]
        else:
            res['%s_time_s' % lab] = np.NaN
            res['%s' % lab] = np.NaN            
            
    # Plot Ca
    ax[Ca_row-1].plot(d['time'], d['Ca'], '-')
    ax[Ca_row-1].plot(res['stim_time_s'], res['Ca_initial'], 'ko',
                      fillstyle='none')
    ax[Ca_row-1].plot(res['Ca_max_time_s'], res['Ca_max'], 'ko',
                      fillstyle='none')
    # Plot times    
    for t in thresh:
        y_lab = ('Ca_up_%.0f' % (100 * t))
        x_lab = '%s_time_s' % y_lab
        ax[Ca_row-1].plot(res[x_lab], res[y_lab], 'rs', fillstyle='none')
        
    for t in thresh:
        y_lab = ('Ca_down_%.0f' % (100 * t))
        x_lab = '%s_time_s' % y_lab
        ax[Ca_row-1].plot(res[x_lab], res[y_lab], 'rs', fillstyle='none')
    
    
    # Plot dCadt
    ax[dCadt_row-1].plot(d['time'], d['dCadt'], '-')
    ax[dCadt_row-1].plot(res['dCadt_max_time_s'], res['dCadt_max'], 'bo',
                         fillstyle='none')
    ax[dCadt_row-1].plot(res['dCadt_min_time_s'], res['dCadt_min'], 'bo',
                         fillstyle='none')
    
    # Plot hs_length
    ax[hsl_row-1].plot(d['time'], d['hs_1_length'], '-')
    ax[hsl_row-1].plot(res['stim_time_s'], res['hs_length_initial'], 'ko',
                      fillstyle='none')
    ax[hsl_row-1].plot(res['hs_length_min_time_s'], res['hs_length_min'], 'ko',
                      fillstyle='none')

    for t in thresh:
        y_lab = ('hs_length_down_%.0f' % (100 * t))
        x_lab = '%s_time_s' % y_lab
        ax[hsl_row-1].plot(res[x_lab], res[y_lab], 'rs', fillstyle='none')
        
    for t in thresh:
        y_lab = ('hs_length_up_%.0f' % (100 * t))
        x_lab = '%s_time_s' % y_lab
        ax[hsl_row-1].plot(res[x_lab], res[y_lab], 'rs', fillstyle='none')        

    # Finally plot dhsldt
    ax[d2Cadt2_row - 1].plot(d['time'], d['d2Cadt2'], '-')
    
    # Create an image file string
    output_image_file_string = os.path.join(analysis_folder,
                                       'images',
                                       ('sim_%i.png' % (sim_index + 1)))
    # Make sure the directory exists
    parent_dir = Path(output_image_file_string).parent.absolute()
    if not (os.path.isdir(parent_dir)):
        os.makedirs(parent_dir)
    
    print('Saving analysis_figure to: %s' % output_image_file_string)
    fig.savefig(output_image_file_string, bbox_inches='tight')
    
    # Tidy up
    plt.close(fig)
    
    # Convert the dict to a dataframe and return
    results_df = pd.DataFrame([res])
    
    return results_df
    
def last_point_before(t, y, t_min, t_max, y_thresh, direction):
    
    if (direction == 'below'):
        vi = np.nonzero((t > t_min) &
                        (t < t_max) &
                        (y.to_numpy() < y_thresh))[0]
    else:
        vi = np.nonzero((t > t_min) &
                        (t < t_max) &
                        (y.to_numpy() > y_thresh ))[0]
 
    if (len(vi) > 0):
        vi = vi[-1]
    else:
        vi = -1
    
    return vi

def summary_figure(top_data_folder, sim_file_string, analysis_folder):
    """ Superposes simulations """
    
    # Set up a figure
    no_of_rows = 2
    no_of_cols = 1
        
    Ca_row = 1
    hs_length_row = 2
    
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([3, 7])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
    
    # Get the paths
    data_dirs = [f for f in Path(top_data_folder).iterdir() if f.is_dir()]
    for d in data_dirs:
        dfs = os.path.join(os.path.abspath(str(d)),
                           sim_file_string)
        
        d = pd.read_csv(dfs, sep='\t', na_values=['-nan(ind'])
        
        d['Ca2+'] = np.power(10, -d['hs_1_pCa'])
        
        ax[Ca_row - 1].plot(d['time'], d['Ca2+'], '-')
        ax[Ca_row - 1].set_ylabel('Ca2+\n[M]')
        
        ax[hs_length_row - 1].plot(d['time'], d['hs_1_length'], '-')
        ax[hs_length_row - 1].set_ylabel('Half-sarcomere\nlength\n(µm)')
        ax[hs_length_row - 1].set_xlabel('Time_s')
        
    ofs = os.path.join(analysis_folder, 'summary.png')
    print('Saving summary_figure to: %s' % ofs)
    fig.savefig(ofs, bbox_inches='tight')

    
    
if __name__ == "__main__":

    top_data_folder = '../sim_data'
    sim_file_string = 'sim_output/1/sim_prot_1_r1.txt'

    excel_data_file_string = '../generated/parameter_values.xlsx'
    analysis_folder = '../analysis'

    base_folder = str(Path(sys.argv[0]).parent.absolute())
    analysis_folder = os.path.join(base_folder, analysis_folder)
    excel_data_file_string = os.path.join(base_folder, excel_data_file_string)
    top_data_folder = os.path.join(base_folder, top_data_folder)

    # Clean the analysis folder
    try:
        print('Trying to remove: %s' % analysis_folder)
        shutil.rmtree(analysis_folder, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s' % (analysis_folder, e.strerror))

    if not os.path.isdir(analysis_folder):
        print('Making %s')
        os.makedirs(analysis_folder)

    pair_plot(excel_data_file_string, analysis_folder)
    
    summary_figure(top_data_folder, sim_file_string, analysis_folder)
    
    analyze_simulations(top_data_folder, sim_file_string, analysis_folder)
