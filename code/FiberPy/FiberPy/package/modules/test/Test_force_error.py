# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 20:07:43 2020

@author: Sarah
"""


# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 09:22:04 2020

@author: Sarah
"""

import os
import subprocess
import Force_balance_from_dump as fb
import numpy as np
import matplotlib.pyplot as plt
import json

PLOT_FIGURES = True 

ROOT = os.path.dirname(__file__)
FIBERSIM_ROOT = os.path.join(ROOT, "..", "..","..", "..", "..", "..", "..", "FiberSim2020_test")
FIBERSIM_EXE = os.path.join(FIBERSIM_ROOT, "bin","FiberSim.exe")
INPUT_DIR = os.path.join(ROOT, "input_files")
MODEL_FILE = os.path.join(INPUT_DIR, "model_force.json")
PROTOCOL_FILE = os.path.join(INPUT_DIR, "pCa45_protocol.txt")

OPTION_DICT = {
  "options": {
    "log_folder": None,
    "dump_hs_status": 1,
    "particle_no": 43,
    "no_of_repeats": 5,
    "max_executable_time": 4800,
    "multithreading": 1,
    "dump_mode": 1,
    "delta_G0": -54e-21,
    "max_rate": 1e4,
    "x_pos_rel_tol": 1e-3
  }
}

precision = [4,8,12,15]

# Run a FiberSim simulation for different dump precisions
force_errors = []

for prec in precision: 

    OPTION_DICT["options"]["log_folder"] = f"output_files/dump_{prec}"
    OPTION_DICT["options"]["dump_precision"] = prec
    OPTION_FILE = os.path.join(INPUT_DIR, "options_test.json")
    with open(OPTION_FILE, 'w') as f:
        json.dump(OPTION_DICT, f)
    
    OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp_{prec}")
    
    HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump_{prec}/hs_status")
    ERROR_FOLDER = os.path.join(ROOT, "output_files")


    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)

# Calculate the thin and thick filaments errors, store it in a .txt file
    
    ref_error = 0.0


    for dump_file in os.listdir(HS_STATUS_FOLDER):
        dump_file = os.path.join(HS_STATUS_FOLDER, dump_file)
        #print (f"Validating: {dump_file}")
    
        error_force = fb.Get_error_force(dump_file)
    
        ref_error = max( (ref_error, error_force) )


    force_errors.append(ref_error)

   
error_file_name = os.path.join(ERROR_FOLDER, "force_errors.txt")
error_storage = np.asarray(force_errors)
np.savetxt(error_file_name, error_storage)

# Plot the error as a function of the dump precision
 
if PLOT_FIGURES:
    

    ERROR_FOLDER = os.path.join(ROOT, "output_files")
    
    error_str = "force_errors.txt"
    error_file = os.path.join(ERROR_FOLDER, error_str)

    this_error = np.loadtxt(error_file)


    plt.figure()
    plt.plot(precision, this_error)
    plt.ylabel("Force Error [mN/mm²]")
    plt.xlabel("Dump Precision [number of digits]")
    plt.yscale("log")
    plt.show()
    







