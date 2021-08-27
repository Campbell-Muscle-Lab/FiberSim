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

ROOT = os.path.dirname(__file__)


def generate_models():
    
    iso_prop = [0.01, 0.1, 0.5, 0.8, 0.99]
    
    model_template = os.path.join(ROOT, "model_template.json")
    
    for i, c in enumerate(iso_prop):
    
        with open(model_template, "r") as jsonFile:
            data = json.load(jsonFile)

        data["m_parameters"]["m_isotype_proportions"] = [1-c, c]
        
        model_string = os.path.join('../', 'sim_input','%i' % i, 'model_%i.json' %i)
        
        # Check the directory exists
        dir_name = os.path.dirname(model_string)
        if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)

        with open(model_string, "w") as jsonFile:
            json.dump(data, jsonFile)
            
        # write output handler as well
        write_output_handler(i)
        
    # Create a job_struct
    job = []
    job_file_string = os.path.join('..', "job_struct.json")

    # Loop through relative force values
    for i, c in enumerate(iso_prop):

        output_handler_file_string = write_output_handler(i)
        model_string = os.path.join('sim_input','%i' % i, 'model_%i.json' %i)
        prot_file_string = os.path.join('sim_input', 'pCa45_protocol.txt')

        # Add in a job
        j=dict()
        j['relative_to'] = 'this_file'
        j['model_file'] = model_string
        j['options_file'] = 'sim_input/sim_options.json'
        j['protocol_file'] = prot_file_string
        j['results_file'] = os.path.join('sim_output', '%i' % i, 'results_%i.txt' % i)
                                         
        j['output_handler_file'] = output_handler_file_string

        job.append(j)

    with open(job_file_string,'w') as f:
        json.dump(job, f, indent=4)

def write_output_handler(j):
    # Writes an output handler file

    # Generate a file string
    file_string = os.path.join('sim_input', '%i' % j, 'output_handler_%i.json' % j)

    # Create the json struct
    ti = []
    ti_j = dict()
    ti_j['relative_to'] = 'this_file'
    ti_j['template_file_string'] = '../../template/template_summary.json'
    ti_j['output_file_string'] = os.path.join('../../sim_output', '%i' %j, 'summary_%i.png' %j)
                                              
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
        
    return file_string


if __name__ == "__main__":
    generate_models()
