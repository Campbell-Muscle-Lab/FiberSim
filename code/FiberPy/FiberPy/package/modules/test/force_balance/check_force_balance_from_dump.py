# -*- coding: utf-8 -*-
"""
Created on Wed May 26 14:00:59 2021

@author: sako231
"""


import os, sys
import json

from natsort import natsorted

import path_definitions as pd
import force_balance as fb

import matplotlib.pyplot as plt

ROOT = os.path.dirname(__file__)
PACKAGE_ROOT = os.path.realpath(os.path.join(ROOT, "..", "..", "..", ".."))
sys.path.append(PACKAGE_ROOT)

from package.modules.half_sarcomere import half_sarcomere
from package.modules.batch import batch as b

def error_and_dump_precision(precision_array, run_sim = True, output_dir = []):
    
    if run_sim:
        
        for prec in precision_array:
            
            # Set dump_precision value in option file
            
            with open(pd.OPTION_FILE, 'r') as f:
                opt = json.load(f)
                
                opt["options"]["dump_precision"] = prec
                
                opt["options"]["status_files"]["status_folder"] =  f"../sim_output/hs_{prec}"
                
            with open(pd.OPTION_FILE, 'w') as f:
                
                json.dump(opt, f)
                
            # Run FiberSim simulation
            
            b.run_batch(pd.BATCH_FILE)
      
    
    max_thick_node_err = []
    max_thin_node_err = []
    
    total_force_err = []
    
    for prec in precision_array:
        
        dump_folder = pd.HS_STATUS_FOLDER + f"_{prec}"
    
        hs_file = []

        for filenames in os.listdir(dump_folder):
            
            filenames = os.path.join(dump_folder, filenames)
            hs_file.append(filenames)
    
        hs_file = natsorted(hs_file) # sorting the dump files    
    
        hs = half_sarcomere.half_sarcomere(hs_file[-1]) # get last dump_file

        thick_err, thin_err = fb.get_hs_thin_and_thick_errors(hs)
        
        max_thick_node_err.append(max(thick_err))
        
        max_thin_node_err.append(max(thin_err))

        total_force_err.append(fb.check_total_force(hs))
        
    # print(max_thick_node_err)
    # print(max_thin_node_err)
        
    plt.figure()
    plt.plot(precision_array, total_force_err, color = "tab:red")
    plt.xlabel("dump_precision")
    plt.xticks(precision_array)
    plt.title("Total force error (%)")
    plt.yscale("log")
    
    if output_dir:    
        output_file = os.path.join(output_dir, "error_force.png")
        plt.savefig(output_file, dpi = 300)
    
    plt.figure()
    plt.plot(precision_array, max_thick_node_err, color = "tab:red", label = "max thick node error")
    plt.plot(precision_array, max_thin_node_err, color = "tab:blue", label = "max thin node error")
    plt.xlabel("dump precision (# digits)")
    plt.xticks(precision_array)
    plt.yscale("log")
    plt.ylabel("[nm]")
    plt.legend()
    
    if output_dir:  
        output_file = os.path.join(output_dir, "error_filaments.png")
        plt.savefig(output_file, dpi = 300)
        
error_and_dump_precision([4,6,8], run_sim = True, output_dir = [])