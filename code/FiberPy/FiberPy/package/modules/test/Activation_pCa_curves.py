# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 10:01:54 2020

@author: Sarah
"""

import os
import subprocess
import numpy as np
import json

import matplotlib.pyplot as plt


# General path definitions

ROOT = os.path.dirname(__file__)
FIBERSIM_ROOT = os.path.join(ROOT, "..", "..","..", "..", "..", "..", "..", "FiberSim")
FIBERSIM_EXE = os.path.join(FIBERSIM_ROOT, "bin","FiberSim.exe")
INPUT_DIR = os.path.join(ROOT, "input_files")
MODEL_FILE = os.path.join(INPUT_DIR, "model_thin_kinetics.json")
OPTION_FILE = os.path.join(INPUT_DIR, "options.json")
OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp")
RESULT_FILE = os.path.join(OUTPUT_DIR, "results.txt")


def pCa_list():
    
    pCa_list = [40]
    
    for i in range(0,21):
        
        pCa_list.append(45 + i)
        
    pCa_list.append(70)
    
    #print(pCa_list)
        
    return pCa_list
    
    
def create_protocol_files(t_step, t_max, pCa_list):
    
    n_t_step = int(t_max/(t_step))
                   
    dt = np.zeros(n_t_step)
    calcium = np.zeros(n_t_step)
    dhsl = np.zeros(n_t_step)
    mode = np.zeros(n_t_step)
                                      
    for ind, pCa in enumerate(pCa_list):   
        
        filename = f"input_files/activation_pCa_curves/pCa_{pCa}.txt"
        
        for delta_t in range(0,n_t_step): 
            
            dt[delta_t] = t_step
            calcium[delta_t] = pCa/10
            dhsl[delta_t] = 0.0
            mode[delta_t] = -2.0
            
            with open(filename, 'wb') as fout:
                NEWLINE_SIZE_IN_BYTES = 2 
                np.savetxt(filename, list(zip(dt,calcium,dhsl,mode)), fmt = "%.4f %.1f %.1f %.1f", header = "dt     pCa dhsl mode", comments='')
                fout.seek(0, os.SEEK_END) # Go to the end of the file.
                # Go backwards one byte from the end of the file.
                fout.seek(fout.tell() - NEWLINE_SIZE_IN_BYTES, os.SEEK_SET)
                fout.truncate() # Truncate the file to this point.
                
def fit_pCa_data(x,y):
    """ Fits Hill-curve to x-y data """
    
    from scipy.optimize import curve_fit

    def y_pCa(x_data, pCa_50, n_H, y_min, y_amp):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = y_min + \
                y_amp * (np.power(np.power(10, -x), n_H) /
                (np.power(np.power(10, -x), n_H) + 
                     np.power(np.power(10, -pCa_50), n_H)))
        return y

    popt, pcov = curve_fit(y_pCa, x, y, [6.0, 2, np.amin(y), np.amax(y)])
    
    d = dict()
    d['pCa_50'] = popt[0]
    d['n_H']    = popt[1]
    d['y_min']  = popt[2]
    d['y_max']  = popt[3]
    d['x_fit']  = np.linspace(9.0, 4, 1000)
    d['y_fit']  = y_pCa(d['x_fit'], *popt)

    return d, pcov
    

time_step = 0.001
time_max = 0.05
pCa_list = pCa_list()
create_protocol_files(time_step,time_max,pCa_list)  


calcium = []
activated_bs = []


for pCa in pCa_list:
    
    with open(MODEL_FILE, "r") as jsonFile:
        data = json.load(jsonFile)

        a_k_coop = data["thin_parameters"]["a_k_coop"]
        data["thin_parameters"]["a_k_coop"] = 10
        
        a_k_off = data["thin_parameters"]["a_k_off"]
        data["thin_parameters"]["a_k_off"] = 100
        
        a_k_on = data["thin_parameters"]["a_k_off"]
        data["thin_parameters"]["a_k_on"] = 1e+7
        
    with open(MODEL_FILE, "w") as jsonFile:
        json.dump(data, jsonFile) 
    
    PROTOCOL_FILE = os.path.join(INPUT_DIR, f"activation_pCa_curves/pCa_{pCa}.txt")
    
    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)
    
    result = np.loadtxt(RESULT_FILE, skiprows=1)   
    calcium_conc = pow(10,-result[0][1])
    calcium.append(calcium_conc)
    active_bs_pop = result[-1][7]
    activated_bs.append(active_bs_pop)
    

filename = f"output_files/activation_pCa.txt"   
np.savetxt(filename,(calcium,activated_bs))
results = np.loadtxt(f"output_files/activation_pCa.txt")

plt.figure()
plt.semilogx(results[0], results[1],'o-')

p_Ca = -np.log10(results[0])
res, pcov = fit_pCa_data(p_Ca,results[1])
x_data = np.power(10,-res["x_fit"])
plt.plot(x_data, res["y_fit"])

n_H = res["n_H"]
n_H_err = np.sqrt(np.diag(pcov))[1]
pCa_50 = res["pCa_50"]
pCa_50_err = np.sqrt(np.diag(pcov))[0]

print("Hill coeffficient = ", res["n_H"])
print("pCa50 = ", res["pCa_50"])


    


