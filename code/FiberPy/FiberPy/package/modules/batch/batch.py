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


def run_batch(json_batch_file_string=[],
              batch_structure=[]):
    """Runs >=1 simulation using multithreading"""

    # We need a batch structure. Load it from file
    # if it is not provided
    if (batch_structure == []):
        print('Running batch from %s' % batch_structure)
        with open(json_batch_file_string, 'r') as f:
            json_data = json.load(f)
            batch_structure = json_data['FiberSim_batch']

    # Pull off the exe path
    exe_structure = batch_structure['FiberSim_exe']
    exe_string = exe_structure['exe_path']
    if (not exe_structure['relative_to']):
        exe_string = os.path.abspath(exe_string)
    elif (exe_structure['relative_to'] == 'this_file'):
        base_directory = Path(json_batch_file_string).parent.absolute()
        exe_string = os.path.join(base_directory, exe_string)
    else:
        base_directory = exe_structure['relative_to']
        exe_string = os.path.join(base_directory, exe_string)

    print('exe_string: %s' % exe_string)

    # Parse the job data into a list of command strings
    job_data = batch_structure['job']
    command_strings = []
    for i,j in enumerate(job_data):
        # Build up command string
        com_string = exe_string
        for f in ['model_file_string', 'options_file_string',
                  'protocol_file_string', 'output_folder']:
            fs = j[f]
            if (not j['relative_to']):
                fs = os.path.abspath(fs)
            elif (j['relative_to'] == 'this_file'):
                base_directory = Path(json_batch_file_string).parent.absolute()
                fs = os.path.join(base_directory, fs)
            else:
                base_directory = j['relative_to']
                fs = os.path.join(base_directory, fs)

            com_string = '%s %s' % (com_string, fs)

        command_strings.append(com_string)
        
        print(command_strings)

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
                if not thread.isAlive():
                    threads.remove(thread)
    
    # At this point we have run all the simulations
# Cycle through the jobs checking for output handlers
                    


def worker(cmd):
    subprocess.call(cmd)


def run_my_batch(json_batch_file_string):
    """Runs FiberSim simulations described in a json batch file"""
    
    print("running batch: %s" % json_batch_file_string)
    
    command_strings = []
    exe_string = "..\\FiberSim.exe"

    with open(json_batch_file_string) as json_file:
        json_data = json.load(json_file)
        FiberSim_batch = json_data['FiberSim_batch']
        job_data = FiberSim_batch['job']
        for i,j in enumerate(job_data):
            model_file_string = j['model_file_string']
            options_file_string = j['options_file_string']
            protocol_file_string = j['protocol_file_string']
            results_file_string = j['results_file_string']

            # Generate a command
            com_string = "%s %s %s %s %s" % (exe_string,
                                             model_file_string,
                                             options_file_string,
                                             protocol_file_string,
                                             results_file_string)
            
            #print(com_string)
            
            command_strings.append(com_string)

    num_processes = multiprocessing.cpu_count()
    my_list = command_strings
    
    threads=[]
    
    while threads or my_list:
        if (len(threads) < num_processes) and my_list:
            t = threading.Thread(target=worker, args=[my_list.pop()])
            t.setDaemon(True)
            t.start()
            threads.append(t)
        else:
            for thread in threads:
                if not thread.isAlive():
                    threads.remove(thread)
    


