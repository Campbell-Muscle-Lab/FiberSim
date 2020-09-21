import os,sys
import subprocess
import numpy as np
import json

import matplotlib.pyplot as plt

from natsort import natsorted

from pathlib import Path


ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

# General path definitions

ROOT = os.path.dirname(__file__)
FIBERSIM_ROOT = os.path.join(ROOT, "..", "..","..", "..", "..", "..", "..", "FiberSim")
FIBERSIM_EXE = os.path.join(FIBERSIM_ROOT, "bin","FiberSim.exe")
INPUT_DIR = os.path.join(ROOT, "input_files")
MODEL_FILE = os.path.join(INPUT_DIR, "model_thin_kinetics.json")
OPTION_FILE = os.path.join(INPUT_DIR, "options.json")
PROTOCOL_FILE = os.path.join(INPUT_DIR, "pCa5_5_protocol.txt")
OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp")
HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump/hs_status")
ERROR_FOLDER = os.path.join(ROOT, "output_files")


def Get_kinetics_from_mod_file(model_file):
    """
    Extract the thin and thick kinetics parameters from the model file and stores them in a dictionary 

    Parameters
    ----------
    model_file : char
        Path to the model file.

    Returns
    -------
    thin_params : dict
        Dictionary where thin kinetics parameters are stored

    """

    # Extract the thin kinetics from the json model file 

    thin_params = {}
    with open(model_file, 'r') as f:
        mod = json.load(f)
    
    thin_params["a_k_on"] = mod["thin_parameters"]["a_k_on"]
    thin_params["a_k_off"] = mod["thin_parameters"]["a_k_off"]
    thin_params["a_k_coop"] = mod["thin_parameters"]["a_k_coop"]    
    thin_params["cb_a_rate"] = mod["m_kinetics"]["scheme"][0]["transition"][0]["rate_parameters"][0]
    
    
    return thin_params      

def off_rate_bs_no_CB():
    """Check that RUs with bound CBs are unable to deactivate"""
    
    #Run a FiberSim simulation 
    
    MODEL_FILE = os.path.join(INPUT_DIR, "model_thin_kinetics_no_cb.json")
    PROTOCOL_FILE = os.path.join(INPUT_DIR, "varying_pCa_protocol.txt")

    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)

    hs_file = []

    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # naturally sorting the dump_files

    # Extract the time step, the number of time steps and calcium concentration from the protocol file

    protocol = np.loadtxt(PROTOCOL_FILE, skiprows=1)
    num_time_step = protocol.shape[0]
    
    # Initialize the arrays for storing the numbers of transitions

    trans_off = 0
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])

    for ind in range(1,num_time_step-1):
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
    
        for i,thin_fil_0 in enumerate(hs_0.thin_fil):
            thin_fil_1 = hs_1.thin_fil[i]
        
            for unit_ind, unit_state in enumerate(thin_fil_0.unit_status):
                
                unit_state_0 = thin_fil_0.unit_status[unit_ind]
                unit_state_1 = thin_fil_1.unit_status[unit_ind]                
                
                bs_idx_start = (unit_ind // 2)*7*2 + (unit_ind %2) 
                
                unit_occupied = 0
                              
                if unit_state_0 == 1: # Unit is activated                   
                    for j in range(0,7): # Check for bound CBs
                        bs_idx = bs_idx_start + 2*j
                        if thin_fil_0.bound_to_m_f[bs_idx] != -1 : #or thin_fil_1.bound_to_m_f[bs_idx] != -1:
                            unit_occupied = 1
                            break
                    if unit_occupied == 1:
                        if unit_state_1 == 0: # Unit deactivates despise the presence of bound CB
                            trans_off += 1                       
                            
        hs_0 = hs_1                       
                  
    return trans_off

def off_rate_bs_with_CB():
    """Check the off-rate of RUs without bound CBs"""
    
    #Run a FiberSim simulation 
    
    MODEL_FILE = os.path.join(INPUT_DIR, "model_thin_kinetics_cb.json")
    PROTOCOL_FILE = os.path.join(INPUT_DIR, "varying_pCa_protocol.txt")

    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)

    hs_file = []

    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # naturally sorting the dump_files

    # Extract the time step, the number of time steps and calcium concentration from the protocol file

    protocol = np.loadtxt(PROTOCOL_FILE, skiprows=1)
    num_time_step = protocol.shape[0]
    time_step = protocol[0][0]
    
    # Initialize the arrays for storing the numbers of transitions
    
    pot_off = 0
    trans_off = 0
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])

    for ind in range(1,num_time_step-1):
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
    
        for i,thin_fil_0 in enumerate(hs_0.thin_fil):
            thin_fil_1 = hs_1.thin_fil[i]
        
            for unit_ind, unit_state in enumerate(thin_fil_0.unit_status):
                
                unit_state_0 = thin_fil_0.unit_status[unit_ind]
                unit_state_1 = thin_fil_1.unit_status[unit_ind]                
                
                bs_idx_start = (unit_ind // 2)*7*2 + (unit_ind %2) 
                
                unit_occupied = 0
                              
                if unit_state_0 == 1: # Unit is activated                   
                    for j in range(0,7): # Check for bound CBs
                        bs_idx = bs_idx_start + 2*j
                        if thin_fil_0.bound_to_m_f[bs_idx] != -1 : #or thin_fil_1.bound_to_m_f[bs_idx] != -1:
                            unit_occupied = 1
                            break
                    if unit_occupied == 0:
                        pot_off += 1
                        if unit_state_1 == 0: # Unit deactivates 
                            trans_off += 1                       
                            
        hs_0 = hs_1

            
    prob_off = trans_off/pot_off
    rate_off = -np.log(1.0 - prob_off) / time_step 
    
    z = 1.96
    
    N_off = pot_off + z**2
    p_off = (trans_off + z**2/2)/N_off
    
    conf_inter_off = p_off + 2*z*np.sqrt(p_off * (1 - p_off)/N_off)
    
    rate_conf_inter_off = -np.log(1.0 -conf_inter_off) / time_step
                  
    return rate_off, rate_conf_inter_off

def Compute_a_kinetics_no_coop():
    """Approximates the governing rate constants for thin filament kinetics when cooperativity is OFF"""
    
    #Run a FiberSim simulation 

    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)

    hs_file = []

    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # naturally sorting the dump_files

    # Extract the time step, the number of time steps and calcium concentration from the protocol file

    protocol = np.loadtxt(PROTOCOL_FILE, skiprows=1)
    num_time_step = protocol.shape[0]
    time_step = protocol[0][0]
    pCa = protocol[0][1]
    calcium = pow(10,-pCa)
    
    # Initialize the arrays for storing the numbers of transitions

    trans_on = 0
    trans_off = 0
    pot_on = 0
    pot_off = 0
    
    blocked_units = 0

    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])

    for ind in range(1,num_time_step-1):
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
    
        for i,thin_fil_0 in enumerate(hs_0.thin_fil):
            thin_fil_1 = hs_1.thin_fil[i]
        
            for unit_ind, unit_state in enumerate(thin_fil_0.unit_status):
                
                unit_state_0 = thin_fil_0.unit_status[unit_ind]
                unit_state_1 = thin_fil_1.unit_status[unit_ind]                
                
                bs_idx_start = (unit_ind // 2)*7*2 + (unit_ind %2) 
                              
                if unit_state_0 == 0:                   
                    pot_on += 1
                    if unit_state_1 == 1: # Unit activates
                        trans_on +=1                        
                else:
                    # Unit may deactivate only if there is no CB bound within the unit
                    unit_occupied = 0
                    for j in range(0,7):
                        bs_idx = bs_idx_start + 2*j
                        if thin_fil_0.bound_to_m_f[bs_idx] != -1 or thin_fil_1.bound_to_m_f[bs_idx] != -1:
                            unit_occupied = 1
                            blocked_units += 1
                            break
                    if unit_occupied == 0:
                        pot_off += 1
                        if unit_state_1 == 0: # Unit deactivates
                            trans_off += 1      
                            
        hs_0 = hs_1
        
    prob_on = trans_on/pot_on
            
    prob_off = trans_off/pot_off
    
    print("# trans off = ", trans_off)
    print("# pot off = ", pot_off)
    print("# blocked units = ", blocked_units)
        
       
    rate_on = -np.log(1.0 - prob_on) / time_step / calcium    
    rate_off = -np.log(1.0 - prob_off) / time_step 
    
    z = 1.96
    N_on = pot_on + z**2
    p_on = (trans_on + z**2/2)/N_on
    
    conf_inter_on = p_on + 2*z*np.sqrt(p_on * (1 - p_on)/N_on)
    
    rate_conf_inter_on = -np.log(1.0 -conf_inter_on) / time_step / calcium
    
    N_off = pot_off + z**2
    p_off = (trans_off + z**2/2)/N_off
    
    conf_inter_off = p_off + 2*z*np.sqrt(p_off * (1 - p_off)/N_off)
    
    rate_conf_inter_off = -np.log(1.0 -conf_inter_off) / time_step
                  
    return rate_on, rate_off, rate_conf_inter_on, rate_conf_inter_off


def Compute_a_kinetics_coop():
    """Approximates the governing rate constants for thin filament kinetics when cooperativity is ON """
    
    #Run a FiberSim simulation 

    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)

    hs_file = []

    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # naturally sorting the dump_files

    # Extract the time step, the number of time steps and calcium concentration from the protocol file

    protocol = np.loadtxt(PROTOCOL_FILE, skiprows=1)
    num_time_step = protocol.shape[0]
    time_step = protocol[0][0]
    pCa = protocol[0][1]
    calcium = pow(10,-pCa)
    
    # Initialize the arrays for storing the numbers of transitions

    trans_on = [0,0,0]
    trans_off = [0,0,0]
    pot_on = [0,0,0]
    pot_off = [0,0,0]
    prob_on = [0,0,0]
    prob_off = [0,0,0]
    rate_on = np.zeros(3)
    rate_off = np.zeros(3)
    rate_conf_inter_on = np.zeros(3)
    rate_conf_inter_off = np.zeros(3)
    
    blocked_units = [0,0,0]

    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])
    
    #no_coop_count =[0,0,0]

    for ind in range(1,num_time_step-1):
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
    
        for i,thin_fil_0 in enumerate(hs_0.thin_fil):
            thin_fil_1 = hs_1.thin_fil[i]
        
            for unit_ind, unit_state in enumerate(thin_fil_0.unit_status):
                
                unit_state_0 = thin_fil_0.unit_status[unit_ind]
                unit_state_1 = thin_fil_1.unit_status[unit_ind]                
                
                bs_idx_start = (unit_ind // 2)*7*2 + (unit_ind %2) 
                              
                # Evaluate the neighbouring units status
                
                coop_left = 0
                coop_right = 0
                coop_idx = 0
                
                if unit_ind > 1:
                    coop_left += thin_fil_0.unit_status[unit_ind - 2]                                        
                    
                if unit_ind < 52:
                    coop_right += thin_fil_0.unit_status[unit_ind + 2]  
                                                        
                coop_idx = coop_left + coop_right
                
                #no_coop_count[coop_idx] += 1
                
                       
                if unit_state_0 == 0:                   
                    pot_on[coop_idx] += 1
                    if unit_state_1 == 1: # Unit activates
                        trans_on[coop_idx] +=1 
                else:
                    pot_off[coop_idx] += 1
                    if unit_state_1 == 0:
                        trans_off[coop_idx] += 1
   
                    # Unit may deactivate only if there is no CB bound within the unit
                    unit_occupied = 0
                    for j in range(0,7):
                        bs_idx = bs_idx_start + 2*j                        
                        if thin_fil_0.bound_to_m_f[bs_idx] != -1:
                            blocked_units[coop_idx] += 1
                            unit_occupied = 1
                            break
                    if unit_occupied == 0:
                        pot_off[coop_idx] += 1
                        if unit_state_1 == 0: # Unit deactivates
                            trans_off[coop_idx] += 1  
     
        hs_0 = hs_1
        
    print("# trans off = ", trans_off)
    print("# pot off = ", pot_off)
    print("# blocked units = ", blocked_units)
              
    for ind in range(0,3):


        if trans_off[ind] != 0:        
            prob_off[ind] = trans_off[ind]/pot_off[ind]
        else:
            print(f"potential number of on-to-off transition is null for cooperativity state n°{ind}")
            prob_off[ind] = 0
            
        if trans_on[ind] != 0:        
            prob_on[ind] = trans_on[ind]/pot_on[ind]
        else:
            print(f"potential number of off-to-on transition is null for cooperativity state n°{ind}")
            prob_on[ind] = 0
    
        rate_on[ind] = -np.log(1.0 - prob_on[ind]) / time_step / calcium    
        rate_off[ind] = -np.log(1.0 - prob_off[ind]) / time_step 
        
        z = 1.96
        
        if rate_on[ind] != 0:
      
            N_on = pot_on[ind] + z**2/2
            p_on = (trans_on[ind] + z**2/2)/N_on
            
            conf_inter_on = p_on + 2*z*np.sqrt(p_on * (1 - p_on)/N_on)
            
            rate_conf_inter_on[ind] = -np.log(1.0 -conf_inter_on) / time_step / calcium
        
        if rate_off[ind] != 0:
            
            N_off = pot_off[ind] + z**2
            p_off = (trans_off[ind] + z**2/2)/N_off
        
            conf_inter_off = p_off + 2*z*np.sqrt(p_off * (1 - p_off)/N_off)
    
        rate_conf_inter_off[ind] = -np.log(1.0 -conf_inter_off) / time_step
              
    return rate_on, rate_off, rate_conf_inter_on, rate_conf_inter_off


def thin_kinetics_analysis():

    #if Path('output_files/thin_kinetics_validation.txt').is_file() == False:
    #    filename = "output_files/thin_kinetics_validation.txt"
    
    param = Get_kinetics_from_mod_file(MODEL_FILE)
    
    protocol = np.loadtxt(PROTOCOL_FILE, skiprows=1)
    num_time_step = protocol.shape[0]
    time_step = protocol[0][0]
    pCa = protocol[0][1]
    
    CB_state = "CB OFF"
    
    if param["cb_a_rate"] != 0:
        
        CB_state = "CB ON"
    
    
    if param["a_k_coop"] != 0:
            
        print("Cooperativity implemented")
        
        a,b,c,d = Compute_a_kinetics_coop()
    
        theo_on = [param['a_k_on'], param['a_k_on']*(1+param['a_k_coop']), param['a_k_on']*(1+2*param['a_k_coop'])]
        theo_off = [param['a_k_off'], param['a_k_off']*(1+param['a_k_coop']), param['a_k_off']*(1+2*param['a_k_coop'])]
        
        #filename = f"output_files/thin_kinetics_validation_{pCa}.txt"    
        #np.savetxt(filename, (a,b,c,d,theo_on, theo_off))
        
        
       
        print(f"Activation rate (no coop) = {a[0]:,.2f} (% of model value : {100*a[0]/param['a_k_on']:.2f}) \n"
          f"Activation rate (coop_1) = {a[1]:,.2f} (% of model value : {100*a[1]/(param['a_k_on']*(1+param['a_k_coop'])):.2f}) \n"
          f"Activation rate (coop_2) = {a[2]:,.2f} (% of model value : {100*a[2]/(param['a_k_on']*(1+2*param['a_k_coop'])):.2f}) \n")
    
        print(f"Deactivation rate (no coop) = {b[2]:,.2f} (% of model value : {100*b[2]/param['a_k_off']:.2f}) \n"
          f"Deactivation rate (coop_1) = {b[1]:,.2f} (% of model value : {100*b[1]/(param['a_k_off']*(1+param['a_k_coop'])):.2f}) \n"
          f"Deactivation rate (coop_2) = {b[0]:,.2f} (% of model value : {100*b[0]/(param['a_k_off']*(1+2*param['a_k_coop'])):.2f}) \n")
    
        
        x= [0,1,2]
        
        fig, ax = plt.subplots()
        ax.bar(x, a,
           yerr=c/2,
           color=['#35838d', '#d64d55', '#960f16'],
           align='center',
           alpha=0.5,
           ecolor='black',
           capsize=10)
        labels = ["No Coop","Coop_1","Coop_2"]
        ax.plot([-0.5, 0.5], [theo_on[0], theo_on[0]], "--", color = '#35838d')
        ax.plot([0.5, 1.5], [theo_on[1], theo_on[1]], "--", color = '#d64d55' )
        ax.plot([1.5,2.5], [theo_on[2], theo_on[2]], "--", color = '#960f16' )
        ax.set_xticks([0,1,2])
        ax.set_title(f"Activation rate (k_coop = {param['a_k_coop']}) - {CB_state} - pCa = {pCa}")
        ax.set_yticks([])
        ax.set_xticklabels(labels)
        plt.show()
        
        fig, ax = plt.subplots()
        ax.bar(x, b,
           yerr=d/2,
           color=['#960f16','#d64d55','#35838d'],
           align='center',
           alpha=0.5,
           ecolor='black',
           capsize=10)
        labels = ["Coop_2","Coop_1","No Coop"]
        ax.plot([-0.5, 0.5], [theo_off[2], theo_off[2]], "--", color = '#960f16')
        ax.plot([0.5, 1.5], [theo_off[1], theo_off[1]], "--", color = '#d64d55')
        ax.plot([1.5,2.5], [theo_off[0], theo_off[0]], "--", color = '#35838d')
        ax.set_xticks([0,1,2])
        ax.set_title(f"Deactivation rate (k_coop = {param['a_k_coop']}) - {CB_state} - pCa = {pCa}")
        ax.set_yticks([])
        ax.set_xticklabels(labels)
        plt.show()
             
    else:
        
        print("Cooperativity not implemented")
        
        a,b,c,d = Compute_a_kinetics_no_coop()
    
        print(f"Calculated Activation rate = {a:,.2f} (% of model value : {100*a/param['a_k_on']:.2f})")
    
        print(f"Calculated Deactivation rate = {b:.2f} (% of model value : {100*b/param['a_k_off']:.2f})")
        
        
        x= [1,2]    
        rate = [a,b]
        err = [c,d]
        
        print(err)
        
        plt.figure()    
        plt.bar(x[0],rate[0], yerr = err[0]/2,color='#960f16')
        plt.xticks([1],["Activation"])
        plt.plot([0.5, 1.5], [param['a_k_on'], param['a_k_on']], "--", color = '#960f16')
        plt.xlim(0,2)
        plt.yticks([])
        plt.title(f"No cooperativity - {CB_state} - pCa = {pCa}")
        
        plt.figure()    
        plt.bar(x[1],rate[1], yerr = [[err[1]/2],[err[1]/2]], color='#35838d')
        plt.xticks([2],["Deactivation"])
        plt.plot([1.5, 2.5], [param['a_k_off'], param['a_k_off']], "--", color = '#35838d')
        plt.xlim(1,3)
        plt.yticks([])
        plt.title(f"No cooperativity - {CB_state} - pCa = {pCa}")
        plt.show()
        
def check_CB_impact_on_rate():
    
    x = [0,1]
        
    rate_off_no_cb = off_rate_bs_no_CB()
    
    rate_off_cb, errors = off_rate_bs_with_CB()
    
    fig, ax = plt.subplots()
      
    ax.bar(x[1],rate_off_cb, yerr = errors/2,color='#960f16', align='center',alpha=0.5,ecolor='black',capsize=10)
    ax.plot([0, 2], [rate_off_no_cb, rate_off_no_cb], color = '#35838d', label = "RUs with bound CBs")
    ax.plot([0.5, 1.5], [100, 100], "--", color = '#960f16', label = "RUs without bound CBs")
    ax.set_xlim(0,2)
    ax.set_xticks([])
    ax.set_yticks([0,100])
    ax.set_title(f"Deactivation rate - varying pCa")
    ax.legend()
        

check_CB_impact_on_rate()
#thin_kinetics_analysis()
 

