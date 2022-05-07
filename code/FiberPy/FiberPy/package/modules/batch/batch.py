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

from ..analysis import atp_cons


def run_batch(json_batch_file_string=[],
              batch_structure=[],
              figures_only = False):
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
    
    # Check the batch to see if max threads have been specified
    if ('max_threads' in batch_structure):
        requested_max_threads = batch_structure['max_threads']
    else:
        requested_max_threads = float("inf")
        
    # Get max threads available
    available_threads = multiprocessing.cpu_count()-1

    # Set processes to mininmum of requested and available
    num_processes = int(min([requested_max_threads, available_threads]))
    print('Running batch using %i threads' % num_processes)

    if (not figures_only):
        # Now run the batch
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
        if ('rates' in batch_figures):
            print('Now generating rates figure')
            for fig_data in batch_figures['rates']:
                analyses.create_rates_figure(fig_data,
                                          json_batch_file_string)

        if ('superposed_traces' in batch_figures):
            print('Now generating superposed_traces figure')
            for fig_data in batch_figures['superposed_traces']:
                analyses.create_superposed_traces_figure(
                    fig_data,
                    json_batch_file_string)

        if ('pCa_curves' in batch_figures):
            print('Now generating tension-pCa curves')
            for fig_data in batch_figures['pCa_curves']:
                analyses.create_y_pCa_figure(fig_data,
                                             json_batch_file_string)

        if ('force_velocity' in batch_figures):
            print('Now generating force-velocity curves')
            for fig_data in batch_figures['force_velocity']:
                analyses.create_fv_and_power_figure(fig_data,
                                          json_batch_file_string)

        if ('k_tr_analysis' in batch_figures):
            print('Now generating k_tr_analysis figure')
            for fig_data in batch_figures['k_tr_analysis']:
                analyses.create_k_tr_analysis_figure(fig_data,
                                          json_batch_file_string)

        if ('ktr' in batch_figures):
            print('Now generating ktr curves')
            for fig_data in batch_figures['ktr']:
                analyses.create_ktr_figure(fig_data,
                                          json_batch_file_string)

        if ('superpose_ktr_plots' in batch_figures):
            print('Now generating superpose ktr plots')
            for fig_data in batch_figures['superpose_ktr_plots']:
                analyses.superpose_ktr_plots(fig_data,
                                          json_batch_file_string)

        if ('myotrope' in batch_figures):
            print('Now generating dose response curve')
            for fig_data in batch_figures['myotrope']:
                analyses.dose_response(fig_data,
                                          json_batch_file_string)

    # Now see if we have to do a kinetics check
    if ('batch_validation' in batch_structure):
        for validation_data in batch_structure['batch_validation']:
                validation.run_validation(validation_data,json_batch_file_string)

    # Now see if we have to calculate ATP consumption rate
    if ('ATP_consumption' in batch_structure):
        print('Now calculating ATP consumption')
        for data in batch_structure['ATP_consumption']:
                atp_cons.get_ATP_cons(data,json_batch_file_string)

    print('FiberPy: run_batch() closing correctly')

def worker(cmd):
    subprocess.call(cmd)

def run_multiple_batch(json_multiple_batch_file_string):
    """Runs multiple batch_files """

    # Load the multiple batches structure

    print('Loading multiple batches from: %s' % json_multiple_batch_file_string)
    with open(json_multiple_batch_file_string, 'r') as f:
        json_data = json.load(f)
        batch_list = json_data['FiberSim_multiple_batch']['batch_list']

    # Run every batch from the batch list

    for batch_file in batch_list:
        run_batch(batch_file)
