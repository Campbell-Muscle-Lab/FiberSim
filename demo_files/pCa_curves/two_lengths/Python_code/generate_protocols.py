# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:47:07 2021

@author: kscamp3
"""

import os
import json

import numpy as np
import pandas as pd

from pathlib import Path

def generate_protocols():
    """ Create protocols """
    
    # Variables
    job_file_string = '../job_struct.json'

    pCa_values = [8.0, 6.0, 5.8, 5.7, 5.6, 5.4, 5.2, 4.8]
    
    prot = dict()
    prot['dt'] = 0.0001
    prot['n_points'] = 200
    prot['initial_pCa'] = 8.0
    prot['n_step_pCa'] = 10
    prot['step_pCa'] = 8.0

    # Create a job_struct
    job = []
    
    # Loop through struct
    for sl in range(0, 2):
        for pCa in pCa_values:
    
            prot['curve'] = (sl+1)
            prot['step_pCa'] = pCa
            prot_file_string = write_protocol(prot)
            output_handler_file_string = write_output_handler(prot)
    
            # Add in a job
            j=dict()
            j['relative_to'] = 'this_file'
            j['model_file'] = os.path.join('sim_input',
                                           ('%.0f' % prot['curve']),
                                           ('model_%.0f.json' % prot['curve']))
            j['options_file'] = 'sim_input/sim_options.json'
            j['protocol_file'] = prot_file_string
            j['results_file'] = os.path.join('sim_output',
                                             ('%.0f' % prot['curve']),
                                             ('pCa_%.0f_results.txt') % (100 * prot['step_pCa']))
            j['output_handler_file'] = output_handler_file_string

            job.append(j)
            
            print(prot_file_string)

    # Write job struct to file
    b = dict()
    b['FiberSim_batch'] = dict()
    b['FiberSim_batch']['job'] = job

    with open(job_file_string,'w') as f:
        json.dump(b, f, indent=4)


def write_protocol(prot):
    # Write a protocol

    # Generate a file string
    file_string = os.path.join('sim_input',
                               ('%i' % prot['curve']),
                               ('pCa_%.0f' % (100 * prot['step_pCa'])),
                               ('prot_pCa_%.0f.txt') % (100 * prot['step_pCa']))
    
    print(file_string)
    print('ken')
    
    # Now generate a protocol
    dt = prot['dt'] * np.ones(prot['n_points'])
    pCa = prot['initial_pCa'] * np.ones(prot['n_points'])
    pCa[prot['n_step_pCa']::] = prot['step_pCa']
    dhsl = np.zeros(prot['n_points'])
    mode = -2 * np.ones(prot['n_points'])
    
    
    d = {'dt': dt, 'pCa': pCa, 'dhsl': dhsl, 'mode': mode}
    df = pd.DataFrame(data=d)

    # Check the directory exists
    dir_name = os.path.dirname('../' + file_string)
    if (not os.path.isdir(dir_name)):
        os.makedirs(dir_name)

    # Write protocol
    df.to_csv('../' + file_string, sep='\t', index=False)
    
    # Return
    return file_string
    

def write_output_handler(prot):
    # Writes an output handler file

    # Generate a file string
    file_string = os.path.join('sim_input',
                              ('%i' % prot['curve']),
                              ('pCa_%.0f' % (100 * prot['step_pCa'])),
                              ('output_handler_%.0f.json') %
                                  (100 * prot['step_pCa']))

    # Create the json struct
    ti = []
    ti_j = dict()
    ti_j['relative_to'] = 'this_file'
    ti_j['template_file_string'] = '../../../template/template_summary.json'
    ti_j['output_file_string'] = os.path.join('../../../sim_output',
                                              ('%i' % prot['curve']),
                                              ('summary_%.0f.png') % (100 * prot['step_pCa']))
    ti.append(ti_j)

    t = dict()
    t['templated_images'] = ti
    
    print(t)

    # Write struct
    
    # Check the directory exists
    dir_name = os.path.dirname('../' + file_string)
    if (not os.path.isdir(dir_name)):
        os.makedirs(dir_name)

    with open('../' + file_string,'w') as f:
        json.dump(t, f, indent=4)

    print(file_string)

    # Return
    return file_string


if __name__ == "__main__":
    generate_protocols()
