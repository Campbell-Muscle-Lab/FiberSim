# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:20:25 2020

@author: kscamp3
"""

import os
import json
import re

import numpy as np
import pandas as pd

import natsort

from pathlib import Path

def fit_pCa_data(x,y):
    """ Fits Hill-curve to x-y data """
    
    from scipy.optimize import curve_fit

    def y_pCa(x_data, pCa_50, n_H, y_min, y_amp):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = y_min + \
                y_amp * (np.power(np.power(10, -x), n_H) /
                (np.power(np.power(10, -x), n_H) + 
                     np.power(np.power(10, -pCa_50), n_H)))
        return y

    popt, pcov = curve_fit(y_pCa, x, y, [6.0, 2, np.amin(y), np.amax(y)])
    
    d = dict()
    d['pCa_50'] = popt[0]
    d['n_H']    = popt[1]
    d['y_min']  = popt[2]
    d['y_max']  = popt[3]
    d['x_fit']  = np.linspace(9.0, 4, 1000)
    d['y_fit']  = y_pCa(d['x_fit'], *popt)

    return d

def save_figure_to_file(f, im_file_string, dpi=250, verbose=1):
    """ Writes an image to file """
    
    dir_path = os.path.dirname(im_file_string)
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
        
    if (verbose):
        print('Saving figure to: %s' % im_file_string)
        
    f.savefig(im_file_string, dpi=dpi)
    
def collate_sim_data(sim_output_directory,
                     no_of_a_states = 2, no_of_m_states = 3,
                     no_of_c_states = 3):
    """ Collates population data from summary output files """
    
    # Get list of directories
    curve_paths = [x for x in Path(sim_output_directory).iterdir() if x.is_dir()]
    curve_folders = []
    for c in curve_paths:
        curve_folders.append(str(c))

    # Remember to sort the folders into a natural order    
    curve_folders = natsort.natsorted(curve_folders)
    
    # Create the holders
    sim_data = dict()
    sim_data['curve'] = []
    sim_data['pCa'] = []
    sim_data['rep'] = []
    for i in range(no_of_a_states):
        sim_data['a_pop_%i' % i] = []
    for i in range(no_of_m_states):
        sim_data['m_pop_%i' % i] = []
    for i in range(no_of_c_states):
        sim_data['c_pop_%i' % i] = []        

    # Loop through the curve folders, finding results files
    for (ci, cf) in enumerate(curve_folders):
        for file in os.listdir(cf):
            if (file.endswith('.txt')):
                dfs = os.path.join(cf, file)
                print(dfs)
                d = pd.read_csv(dfs, sep='\t')
                
                sim_data['curve'].append((ci+1))
                sim_data['pCa'].append(d['pCa'].iloc[-1])
                sim_data['rep'].append(re.findall(r'\d+', dfs)[-1])
                for i in range(3):
                    if (i==0):
                        field_string = 'a_pop_%i'
                        n = no_of_a_states
                    if (i==1):
                        field_string = 'm_pop_%i'
                        n = no_of_m_states
                    if (i==2):
                        field_string = 'c_pop_%i'
                        n = no_of_c_states
                    for j in range(n):
                        sim_data[field_string % j].append(d[field_string % j].iloc[-1])
                
    collated_data = pd.DataFrame(data=sim_data)
    
    return collated_data
                
                
                
        
