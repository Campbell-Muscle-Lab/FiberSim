# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:23:48 2020

@author: kscamp3
"""

import os
import json
import multiprocessing
import threading
import subprocess

from pathlib import Path

from ..output_handler import output_handler as oh

from ..display import analyses

from ..validation import validation


def run_batch(json_batch_file_string=[],
              batch_structure=[]):
    """Runs >=1 simulation using multithreading"""

    print('FiberPy: run_batch() starting')

    # We need a batch structure. Load it from file
    # if it is not provided
    if (batch_structure == []):
        print('Loading FiberSim batch from: %s' % json_batch_file_string)
        with open(json_batch_file_string, 'r') as f:
            json_data = json.load(f)
            batch_structure = json_data['FiberSim_batch']

    # Pull off the exe path
    exe_structure = batch_structure['FiberCpp_exe']
    exe_string = exe_structure['exe_file']
    if (not exe_structure['relative_to']):
        exe_string = os.path.abspath(exe_string)
    elif (exe_structure['relative_to'] == 'this_file'):
        base_directory = Path(json_batch_file_string).parent.absolute()
        exe_string = os.path.join(base_directory, exe_string)
    else:
        base_directory = exe_structure['relative_to']
        exe_string = os.path.join(base_directory, exe_string)


    # Parse the job data into a list of command strings
    job_data = batch_structure['job']
    command_strings = []
    results_file_strings = []
    for i, j in enumerate(job_data):
        # Build up command string
        com_string = exe_string
        for f in ['model_file', 'options_file',
                  'protocol_file', 'results_file']:
            fs = j[f]
            if (not j['relative_to']):
                fs = os.path.abspath(fs)
            elif (j['relative_to'] == 'this_file'):
                base_directory = Path(json_batch_file_string).parent.absolute()
                fs = os.path.join(base_directory, fs)
            else:
                base_directory = j['relative_to']
                fs = os.path.join(base_directory, fs)
            # Store the results file in case you need it to make a figure
            # using the output_handler system
            if (f == 'results_file'):
                results_file_strings.append(fs)

            com_string = '%s "%s"' % (com_string, fs)

        command_strings.append(com_string)

    # Now run the batch using all but 1 cpi
    num_processes = (multiprocessing.cpu_count() - 1)
    my_list = command_strings

    threads = []
    while threads or my_list:
        if (len(threads) < num_processes) and my_list:
            t = threading.Thread(target=worker, args=[my_list.pop()])
            t.setDaemon(True)
            t.start()
            threads.append(t)
        else:
            for thread in threads:
                if not thread.is_alive():
                    threads.remove(thread)

    # At this point we have run all the simulations
    # Run the output handlers
    for i, j in enumerate(job_data):
        if ('output_handler_file' in j):
            fs = j['output_handler_file']
            if (not j['relative_to']):
                fs = os.path.abspath(fs)
            elif (j['relative_to'] == 'this_file'):
                base_directory = Path(json_batch_file_string).parent.absolute()
                fs = os.path.join(base_directory, fs)
            else:
                base_directory = j['relative_to']
                fs = os.path.join(base_directory, fs)
            oh.output_handler(fs,
                              sim_results_file_string=results_file_strings[i])

    # Now see if we have to make any figures
    if ('batch_figures' in batch_structure):
        batch_figures = batch_structure['batch_figures']

        # Dive into the structure
        if ('pCa_curves' in batch_figures):
            for fig_data in batch_figures['pCa_curves']:
                analyses.create_y_pCa_figure(fig_data,
                                             json_batch_file_string)

        if ('force_velocity' in batch_figures):
            for fig_data in batch_figures['force_velocity']:
                analyses.create_fv_and_power_figure(fig_data,
                                          json_batch_file_string)

    # Now see if we have to do a kinetics check
    if ('batch_validation' in batch_structure):
        for validation_data in batch_structure['batch_validation']:
                validation.run_validation(validation_data,json_batch_file_string)

    print('FiberPy: run_batch() closing correctly')

def worker(cmd):
    subprocess.call(cmd)
