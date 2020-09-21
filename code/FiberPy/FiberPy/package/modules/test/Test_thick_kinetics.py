# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:00:32 2020

@author: Sarah
"""


import os,sys
import subprocess
import numpy as np
import json

from natsort import natsorted

import matplotlib.pyplot as plt


ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

# General path definitions

ROOT = os.path.dirname(__file__)
FIBERSIM_ROOT = os.path.join(ROOT, "..", "..","..", "..", "..", "..", "..", "FiberSim")
FIBERSIM_EXE = os.path.join(FIBERSIM_ROOT, "bin","FiberSim.exe")
INPUT_DIR = os.path.join(ROOT, "input_files")
MODEL_FILE = os.path.join(INPUT_DIR, "model_thick_kinetics.json")
OPTION_FILE = os.path.join(INPUT_DIR, "options.json")
PROTOCOL_FILE = os.path.join(INPUT_DIR, "pCa4_protocol.txt")
OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp")
HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump/hs_status")
ERROR_FOLDER = os.path.join(ROOT, "output_files")
CLASS_FILE = os.path.join(OUTPUT_DIR, f"class.json")

# General intervals definitions

X_MIN = -10
X_MAX = 10
X_STEP = 0.25
NB_INTER = int((X_MAX - X_MIN)/X_STEP)


def Get_m_kinetics_from_mod_file(model_file):
    """
    Extract the transitions data from the model file and stores them in a dictionary 

    Parameters
    ----------
    model_file : char
        Path to the model file.

    Returns
    -------
    m_kinetics : dict
        Dictionary where transitions types and parameters are stored

    """

    # Extract the myosin kinetics from the json model file 

    m_kinetics = {}
    with open(model_file, 'r') as f:
        mod = json.load(f)
    
    m_kinetics["no_of_states"] = mod["m_kinetics"]["no_of_states"]
    
    # Create temporary dictionnaries to store the transitions properties

    dict_trans = {}
    dict_rate = {}
    dict_idx = {}
    bound_types = []
    
    idx = 0
    
    for state in mod["m_kinetics"]["scheme"]:
        num = state["number"]
        dict_trans[num] = {}
        dict_rate[num] = {}
        dict_idx[num] = {}
        
        if state["type"] == 'A':
            bound_types.append(state['number'])
            
        for trans in state["transition"]:
            dict_trans[num][trans["new_state"]] = trans["rate_type"]
            dict_rate[num][trans["new_state"]] = trans["rate_parameters"]
            dict_idx[num][trans["new_state"]] = idx
            idx += 1 

    # Store the transitions types/rate parameters/indices in m_kinetics
        
    m_kinetics["transition_type"] = dict_trans
    m_kinetics["transition_rate"] = dict_rate
    m_kinetics["transition_idx"] = dict_idx
    
    m_kinetics['bound_types'] = bound_types
    
    m_kinetics["no_of_transitions"] = idx
    
    m_kinetics["k_cb"] = mod["m_parameters"]["m_k_cb"]
    
    CLASS_FILE = os.path.join(OUTPUT_DIR, f"class.json")
    with open(CLASS_FILE, 'w') as f:
        json.dump(m_kinetics, f)
    
    return m_kinetics


def Get_rates_values(m_kinetics, output_dir):
    """    
    Save all the rate transitions in a .txt file
    
    Parameters
    ----------
    m_kinetics : dict
        Contains the transition properties
        
    output_dir : char
        Path to the folder where the .txt file is stored

    Returns
    -------
    None.

    """

    stretch = np.arange(X_MIN,X_MAX, X_STEP)
    
    rate_values = [stretch]
    
    #rate_values = {}

    for state in m_kinetics["transition_type"]:
        for new_state in m_kinetics["transition_type"][state]:
            
            trans_type = m_kinetics["transition_type"][state][new_state]
            trans_param = m_kinetics["transition_rate"][state][new_state]
            #idx = m_kinetics["transition_idx"][state][new_state]
            
            if trans_type == "constant":
                
                #rate_values[idx] = [trans_param[0] for x in stretch]
                rate_trans = [trans_param[0] for x in stretch]
                rate_values.append(rate_trans)
                               
            elif trans_type == "gaussian":
                
                k_cb = m_kinetics["k_cb"]
                rate_trans = [trans_param[0]*np.exp(-0.5*k_cb*np.power(x, 2)/(1e18*1.38e-23*310)) for x in stretch]
                rate_values.append(rate_trans)
                                
            elif trans_type == "poly":
                
                k_cb = m_kinetics["k_cb"]
                rate_trans = [trans_param[0] + trans_param[1]*np.power(x, trans_param[2]) for x in stretch]
                rate_values.append(rate_trans)
                
            else:
                raise RuntimeError(f"Transition of type {trans_type} is unknown")
                
    #for element in range(m_kinetics["no_of_transitions"]):
        #array.append(rate_values[element])
    filename = os.path.join(output_dir, f"rate_equations.txt")       
    np.savetxt(filename, np.transpose(rate_values))
    
    
def Get_stretch_interval(stretch):
    """ Provide the interval number associated with a given stretch value
    
    Parameters
    ----------
    stretch : float
        CB stretch

    Returns
    -------
    no_interval : int
        # interval
    
    """
        
    no_interval = int(np.floor((stretch - X_MIN) / X_STEP))
    
    if no_interval >= NB_INTER:
        #print('interval out of bounds')
        no_interval = NB_INTER-1
        
    if no_interval <= 0:
        #print('interval out of bounds')
        no_interval = 0

    return no_interval
    
    

def Compute_m_kinetics_rate():
    """Approximates the governing rate functions of FiberSim as a function of stretch (delta x)"""

    hs_file = []

    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # naturally sorting the dump_files

    # Extract the time step and the number of time steps from the protocol file

    protocol = np.loadtxt(PROTOCOL_FILE, skiprows=1)
    num_time_step = protocol.shape[0]
    time_step = protocol[0][0]

    # Initialize the arrays for storing the numbers of transitions

    complete_transition = np.zeros((m_kinetics["no_of_transitions"], NB_INTER),dtype=int)
    potential_transition = np.zeros((m_kinetics["no_of_transitions"], NB_INTER),dtype=int)
    prob = np.zeros((m_kinetics["no_of_transitions"], NB_INTER),dtype=float)
    calculated_rate = np.zeros((m_kinetics["no_of_transitions"], NB_INTER),dtype=float)


    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])


    for ind in range(1,num_time_step-1):
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
    
        for i,thick_fil_0 in enumerate(hs_0.thick_fil):
            thick_fil_1 = hs_1.thick_fil[i]
        
            for cb_ind, state in enumerate(thick_fil_0.cb_state):
                cb_state_0 = thick_fil_0.cb_state[cb_ind]
                cb_state_1 = thick_fil_1.cb_state[cb_ind]
                
                # Get the attached CB stretch               
                if cb_state_0 in m_kinetics["bound_types"]: 

                    thin_ind = thick_fil_0.cb_bound_to_a_f[cb_ind]
                    bs_ind = thick_fil_0.cb_bound_to_a_n[cb_ind]
                # Get the "potentiel" CB stretch if attachment had occurred
                else:                                       
                    thin_ind = thick_fil_0.cb_nearest_a_f[cb_ind]
                    bs_ind = thick_fil_0.cb_nearest_a_n[cb_ind]
                    
                thin_fil = hs_0.thin_fil[thin_ind]
                stretch = thick_fil_0.cb_x[cb_ind] - thin_fil.bs_x[bs_ind] 
                
                cb_stretch_0 = Get_stretch_interval(stretch) 
                    
                if cb_state_0 != cb_state_1: # A transition has occurred
                    
                    #print("A transition occurred")
                    idx = m_kinetics["transition_idx"][cb_state_0][cb_state_1]                                 
                    complete_transition[idx,cb_stretch_0] += 1
                    potential_transition[idx,cb_stretch_0] += 1
                
                # Check for all the potential transitions that could have lead cb_0 to cb_1 
                for transition in m_kinetics["transition_idx"][cb_state_0]:
                    
                    pot_trans = False
                    cb_state_1 = transition
                    
                    if cb_state_1 in m_kinetics["bound_types"]:
                        pot_trans = True 
                        
                        # Check if the potential binding site on thin_1 is available 
                        
                        thin_ind_1 = thick_fil_1.cb_nearest_a_f[cb_ind] # find the nearest thin filament where attachment could occur
                        
                        bs_ind_0 = thick_fil_0.cb_nearest_a_n[cb_ind] # bs node index at time t   
                        bs_ind_1 = thick_fil_1.cb_nearest_a_n[cb_ind] # bs node index at time t + time_step
                        
                        thin_fil_0 = hs_0.thin_fil[thin_ind_1] # thin fil at time t 
                        thin_fil_1 = hs_1.thin_fil[thin_ind_1] #thin fil at time t + time_step  
                        
                                                
                        if thin_fil_0.bs_state[bs_ind_0] == 0:
                            pot_trans = False # binding site was "off"
                    
                        if thin_fil_0.bound_to_m_f[bs_ind_0] >= 0:
                            pot_trans = False # binding site was "occupied" 
                            
                        if thin_fil_1.bs_state[bs_ind_1] == 0:
                            pot_trans = False # binding site is "off"                             
                            
                        # if binding site is "on", check if there's already an attached cb that is not cb_ind  
                            
                        if thin_fil_1.bound_to_m_f[bs_ind_1] >= 0 and thin_fil_1.bound_to_m_f[bs_ind_1] != thick_fil_0.thick_id:
                            if thin_fil_1.bound_to_m_n[bs_ind_1] >= 0 and thin_fil_1.bound_to_m_n[bs_ind_1] != cb_ind:
                                pot_trans = False
                                
                        if pot_trans:
                            
                            stretch = thick_fil_0.cb_x[cb_ind] - thin_fil_0.bs_x[bs_ind_1]
                            cb_stretch_1 = Get_stretch_interval(stretch)                                                           
                        
                    else:
                        pot_trans = True # Detachment could have occurred 
                        cb_stretch_1  = cb_stretch_0 
                        
                    if pot_trans:
                        idx_pot = m_kinetics["transition_idx"][cb_state_0][cb_state_1]                                                      
                        potential_transition[idx_pot, cb_stretch_1] += 1
                            
        hs_0 = hs_1
        
    for eqn in range(0, m_kinetics["no_of_transitions"]):
        for ind in range(0,NB_INTER-1):
            if potential_transition[eqn,ind] == 0:          
                prob[eqn,ind] = 0
            else:
                prob[eqn,ind] = complete_transition[eqn,ind]/potential_transition[eqn,ind]
                if prob[eqn,ind] > 1:
                    print("probability greater than 1")
                calculated_rate[eqn, ind] = -np.log(1.0 - prob[eqn,ind]) / time_step
                
    x = np.arange(X_MIN,X_MAX, X_STEP)
    
    rates_file = os.path.join(OUTPUT_DIR, f"rate_equations.txt")  
    rate_values = np.loadtxt(rates_file)
    stretches = rate_values[:, 0]
    
    plt.figure()
    plt.plot(x, calculated_rate[0,:])
    plt.ylabel("Rate (s-1)")
    plt.xlabel("CB stretch [nm]")
    plt.title("Attachment rate")
    plt.plot(stretches, rate_values[:, 1])  
    plt.show()  

    plt.figure()
    plt.plot(x, calculated_rate[1,:])
    plt.ylabel("Rate (s-1)")
    plt.xlabel("CB stretch [nm]")
    plt.title("Detachment rate")    
    plt.plot(stretches, rate_values[:, 2])  
    plt.show()          
              
    return calculated_rate

cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
subprocess.call(cmd)
     
m_kinetics = Get_m_kinetics_from_mod_file(MODEL_FILE)  
Get_rates_values(m_kinetics, OUTPUT_DIR)
calc_rate = Compute_m_kinetics_rate()

#print(calc_rate)




    
    

                
        
   









