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


def run_batch(json_batch_file_string=[],
              batch_structure=[]):
    """Runs >=1 simulation using multithreading"""

    # We need a batch structure. Load it from file
    # if it is not provided
    if (batch_structure == []):
        print('Running batch from %s' % batch_structure)
        with open(json_batch_file_string, 'r') as f:
            batch_structure = json.load(f)

    # Now turn the batch structure into command strings
    # Pull off the exe string
    exe_string = os.path.abspath(batch_structure['FiberSim_exe_path'])

    # Now build up the jobs
    command_strings = []
    job_data = batch_structure['job']
    for i, j in enumerate(job_data):
        if ('input_folder' in j):
            folder = j['input_folder']
        else:
            folder = ''
        com_string = exe_string
        for key in list(j.keys()):
            if (key in ['input_folder', 'target_file_string']):
                continue
            com_string = com_string + ' ' + \
                os.path.abspath(os.path.join(folder, j[key]))

        command_strings.append(com_string)
            
    print(command_strings)

    # Now that we have assembled the commands, run the batch
    # using multithreading if useful
    num_processes = multiprocessing.cpu_count()
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
    


