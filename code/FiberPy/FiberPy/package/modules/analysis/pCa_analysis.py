# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 11:26:02 2024

@author: Campbell
"""

import os
import json

import numpy as np
import pandas as pd

from collections import defaultdict

from io import StringIO

from natsort import natsorted

from package.modules.analysis import curve_fitting as cv
from package.modules.utilities import utilities as ut

def pCa_analysis(fig_data, batch_file_string):
    """ Analyzes simulation records to determine force-pCa parameters """
    
    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)
    
    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']

    # Store curve data in lists
    pCa_values = defaultdict(list)
    y_values = defaultdict(list)

    curve_counter = 1
    keep_going = True

    # Create lists to hold data
    curve_index = []
    hs_force = []
    hs_pCa = []
    hs_length = []
        
    # Keep track of max_y
    max_y = -np.inf
    
    # Set up a dictionary to store the curve data
    curve_data = dict()
    curve_data['curve'] = []

    while (keep_going):
        curve_folder = os.path.join(top_data_folder,
                                    ('%i' % curve_counter))

        if os.path.isdir(curve_folder):
            # Find the results files
            for file in os.listdir(curve_folder):
                if (file.endswith('.txt') and not file.startswith('rates')):
                    data_file_string = \
                        os.path.join(curve_folder, file)
                    d = pd.read_csv(data_file_string, delimiter='\t')
                    pCa_values[curve_counter-1].append(d['hs_1_pCa'].iloc[-1])
                    y = d[fig_data['data_field']].iloc[-50:-1].mean() # take the mean force over last 50 points
                    y_values[curve_counter-1].append(y)
                    if (np.amax(y) > max_y):
                        max_y = np.amax(y)
                    
                    # Store data for subsequent output
                    curve_index.append(curve_counter)
                    hs_force.append(y)
                    hs_pCa.append(d['hs_1_pCa'].iloc[-1])
                    hs_length.append(d['hs_1_length'].iloc[-1])

            # Add in curve
            res=cv.fit_pCa_data(pCa_values[curve_counter-1],
                                y_values[curve_counter-1])
            
            # Store data
            d_parameters = pd.DataFrame({'pCa_50': res['pCa_50'],
                                        'n_H': res['n_H'],
                                        'y_min': res['y_min'],
                                        'y_amp': res['y_amp']},
                                        index=[0])
            d_fits = pd.DataFrame({'x_fit': res['x_fit'],
                                    'y_fit': res['y_fit']});
            
            d_curve = pd.concat([d_parameters, d_fits])
            curve_data['curve'].append(d_curve)
           
            # Loop on to the next folder
            curve_counter = curve_counter + 1
    
        else:
            keep_going = False
            
    # Take lists and create a data frame
    r = pd.DataFrame({'curve': curve_index,
                      'hs_pCa': hs_pCa,
                      'hs_force': hs_force,
                      'hs_length': hs_length})
    
    # Save the data as an excel file if required in the batch file
    if('output_data_file_string' in fig_data):
        if (fig_data['relative_to'] == 'this_file'):
            output_file_string = os.path.join(base_folder,
                                              fig_data['output_data_file_string'])
        else:
            output_file_string = fig_data['output_data_file_string']

        print('Writing tension-pCa data to %s' % output_file_string)
        with pd.ExcelWriter(output_file_string, engine='openpyxl') as writer:
            r.to_excel(writer, sheet_name = 'simulation_data', index=False)
            for (i,c) in enumerate(curve_data['curve']):
                c.to_excel(writer,
                           index=False,
                           sheet_name = 'curve_%i' % (i+1))
