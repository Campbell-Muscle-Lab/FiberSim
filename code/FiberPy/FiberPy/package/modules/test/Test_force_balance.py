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
    
    thick_max_error = 0.0
    thin_max_error = 0.0


    for dump_file in os.listdir(HS_STATUS_FOLDER):
        dump_file = os.path.join(HS_STATUS_FOLDER, dump_file)
        #print (f"Validating: {dump_file}")
    
        error_thin, error_thick = fb.Get_error_from_hs_dump_files(dump_file)
    
        thick_max_error = max( (thick_max_error, error_thin) )
        thin_max_error = max( (thin_max_error, error_thick ) )


    error_file_name = os.path.join(ERROR_FOLDER, "hs_errors_dump_{}.txt".format(prec))
    error_storage = np.asarray([thick_max_error, thin_max_error])
    np.savetxt(error_file_name, error_storage)

# Plot the error as a function of the dump precision

if PLOT_FIGURES:
    
    thick_errors = []
    thin_errors = []

    ERROR_FOLDER = os.path.join(ROOT, "output_files")
    
    for prec in precision:
    
        error_str = "hs_errors_dump_{}.txt".format(prec)
        error_file = os.path.join(ERROR_FOLDER, error_str)

        this_error = np.loadtxt(error_file)
        this_thick_err = this_error[0]
        this_thin_err = this_error[1]

        thick_errors.append(this_thick_err)
        thin_errors.append(this_thin_err)

        print(f"Thick errors: {thick_errors}")
        print(f"Thin errors: {thin_errors}")


    plt.figure()
    plt.plot(precision, thick_errors, label="Thick")
    plt.ylabel("Thick Filament Error [nm]")
    plt.xlabel("Dump Precision [number of digits]")
    plt.yscale("log")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(precision, thin_errors, label="Thin")
    plt.ylabel("Thin Filament Error [nm]")
    plt.xlabel("Dump Precision [number of digits]")
    plt.yscale("log")
    plt.legend()
    plt.show()







