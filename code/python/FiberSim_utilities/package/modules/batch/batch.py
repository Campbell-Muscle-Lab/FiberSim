# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:23:48 2020

@author: kscamp3
"""

from pathlib import Path
import json
import multiprocessing
import threading
import subprocess

def run_batch(json_batch_file_string):
    """Runs >=1 simulation using multithreading"""
    
    print("running batch: %s" % json_batch_file_string)
    
    command_strings=[]
    with open(json_batch_file_string) as json_file:
        json_data = json.load(json_file)
        FiberSim_batch = json_data['FiberSim_batch']
        exe_string = Path.resolve(Path.cwd() / FiberSim_batch['FiberSim_exe_path'])
        job_data = FiberSim_batch['job']
        for i,j in enumerate(job_data):
            input_folder_string = j['input_folder']
            model_file_string = Path.resolve(
                Path.cwd() / input_folder_string / j['model_file_string'])
            options_file_string = Path.resolve(
                Path.cwd() / input_folder_string / j['options_file_string'])
            protocol_file_string = Path.resolve(
                Path.cwd() / input_folder_string / j['protocol_file_string'])
            output_folder = Path.resolve(
                Path.cwd() / j['output_folder'])

            # Generate a command
            com_string = "%s %s %s %s %s" % (exe_string,
                                             model_file_string,
                                             options_file_string,
                                             protocol_file_string,
                                             output_folder)
            
            print(com_string)
            
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
    

def worker(cmd):
    subprocess.call(cmd)

