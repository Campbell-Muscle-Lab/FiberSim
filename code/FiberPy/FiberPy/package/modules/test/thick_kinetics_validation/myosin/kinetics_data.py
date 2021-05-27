# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 17:11:49 2021

@author: sako231
"""

import json
import os
import numpy as np
import pandas 

import path_definitions as pd

# Create stretch and node_force bins

X_MIN = -8
X_MAX = 8
X_STEP = 0.5

F_MIN = -0.80
F_MAX = 0.80
F_STEP = 0.05 # same number of bin for stretch and node force

NB_INTER = int((X_MAX - X_MIN)/X_STEP)

A_MIN = 90
A_MAX = 180
A_STEP = 9

NB_A_INTER = int((A_MAX - A_MIN)/A_STEP)

        
def get_stretch_interval(stretch):
    
    # Get the stretch bins for calculating probabilities 
        
    no_interval = int(np.floor((stretch - X_MIN) / X_STEP))
    
    if no_interval >= NB_INTER:
        #print(f'interval out of bounds for stretch = {stretch}')
        no_interval = NB_INTER-1
        
    if no_interval <= 0:
        #print(f'interval out of bounds for stretch = {stretch}')
        no_interval = 0

    return no_interval


def get_node_force_interval(node):
    
    # Get the node force bins for calculating probabilities 
        
    no_interval = int(np.floor((node - F_MIN) / F_STEP))
    
    if no_interval >= NB_INTER:
        no_interval = NB_INTER-1
        
    if no_interval <= 0:
        no_interval = 0

    return no_interval

def get_alignment_interval(angle):
    
    # Get the alignment factor bins for calculating probabilities 
        
    no_interval = int(np.floor((angle - A_MIN) / A_STEP))
    
    if no_interval >= NB_A_INTER:
        no_interval = NB_A_INTER-1
        
    if no_interval <= 0:
        no_interval = 0

    return no_interval


def get_m_kinetics(model_json_file):
    
    # Extract the myosin kinetics from the json model file 

    with open(model_json_file, 'r') as f:
        mod = json.load(f)
           
    m_kinetics = []
    
    for i, isotype in enumerate(mod["m_kinetics"]): # Get data for each isotype
        
        idx = 0
        
        data_scheme = []
        
        for j, state in enumerate(isotype["scheme"]): # Get kinetic scheme for each state
            
            state_data = {}
                        
            state_data["state_number"] = state["number"]
            state_data["state_type"] = state["type"]
            state_data["transition"] = []  # array for the different transitions
                       
            for k, trans in enumerate(state["transition"]): # Get transition data for each new state
                
                trans_data = {}
                
                trans_data["to"] = trans["new_state"]
                trans_data["index"] = idx
                trans_data["rate_type"] = trans["rate_type"]
                trans_data["rate_parameters"] = trans["rate_parameters"]
                
                if state["type"] == 'D' and isotype["scheme"][trans["new_state"] - 1]["type"] == 'A':
                    trans_data["trans_type"] = 'a'
                elif state["type"] == 'A' and isotype["scheme"][trans["new_state"] - 1]["type"] == 'D':
                    trans_data["trans_type"] = 'd'
                elif state["type"] == 'S' and isotype["scheme"][trans["new_state"] - 1]["type"] == 'D':
                    trans_data["trans_type"] = 'srx'
                else:
                    trans_data["trans_type"] = 'x'
                
                idx += 1
                
                state_data["transition"].append(trans_data)
                
            data_scheme.append(state_data)
        
        m_kinetics.append(data_scheme)
        
    # Save the kinetics structure in a JSON file
    outfile = os.path.join(pd.OUTPUT_DIR, "m_kinetics.json")  
    with open(outfile, 'w') as f:
        json.dump(m_kinetics, f)
        
    return m_kinetics

def get_c_kinetics(model_json_file):
    
    # Extract the mybpc kinetics from the json model file 

    with open(model_json_file, 'r') as f:
        mod = json.load(f)
           
    c_kinetics = []
    
    for i, isotype in enumerate(mod["c_kinetics"]): # Get data for each isotype
        
        idx = 0
        
        data_scheme = []
        
        for j, state in enumerate(isotype["scheme"]): # Get kinetic scheme for each state
            
            state_data = {}
                        
            state_data["state_number"] = state["number"]
            state_data["state_type"] = state["type"]
            state_data["transition"] = []  # array for the different transitions
                       
            for k, trans in enumerate(state["transition"]): # Get transition data for each new state
                
                trans_data = {}
                
                trans_data["to"] = trans["new_state"]
                trans_data["index"] = idx
                trans_data["rate_type"] = trans["rate_type"]
                trans_data["rate_parameters"] = trans["rate_parameters"]
                
                if state["type"] == 'D' and isotype["scheme"][trans["new_state"] - 1]["type"] == 'A':
                    trans_data["trans_type"] = 'a'
                elif state["type"] == 'A' and isotype["scheme"][trans["new_state"] - 1]["type"] == 'D':
                    trans_data["trans_type"] = 'd'
                elif state["type"] == 'S':
                    trans_data["trans_type"] = 'srx'
                else:
                    trans_data["trans_type"] = 'x'

                idx += 1
                
                state_data["transition"].append(trans_data)
                
            data_scheme.append(state_data)
        
        c_kinetics.append(data_scheme)
        
    # Save the kinetics structure in a JSON file
    outfile = os.path.join(pd.OUTPUT_DIR, "c_kinetics.json")  
    with open(outfile, 'w') as f:
        json.dump(c_kinetics, f)
        
    return c_kinetics


def calculate_rate_from_m_kinetics(m_kinetics, model_json_file):
    
    # Extract the rate laws from the kinetic data
       
    stretch = np.arange(X_MIN,X_MAX, X_STEP)    
    node_force = np.arange(F_MIN,F_MAX, F_STEP)
    
    with open(model_json_file, 'r') as f:
        mod = json.load(f)

    k_cb = mod["m_parameters"]["m_k_cb"]

    for i, iso in enumerate(m_kinetics):
        
        rate_data = pandas.DataFrame()
        
        rate_data["stretch"] = stretch
        
        rate_data["node_force"] = node_force
    
        for j, state in enumerate(iso):
            
            for new_state in state["transition"]:
                
                trans_type = new_state["rate_type"]
                trans_param = new_state["rate_parameters"]
                index = new_state["index"]
                
                if trans_type == "constant":
                    
                    rate_trans = [trans_param[0] for x in stretch]
                                   
                elif trans_type == "gaussian":
                    
                    rate_trans = [trans_param[0]*np.exp(-0.5 * k_cb * np.power(x + X_STEP/2, 2)/(1e18 * 1.38e-23*310)) for x in stretch]
                                    
                elif trans_type == "poly":
                    
                    rate_trans = [trans_param[0] + trans_param[1] * np.power(x, trans_param[2]) for x in stretch]
                    
                elif trans_type == "force_dependent":
                    
                    rate_trans = [trans_param[0] * (1 + trans_param[1] * max(nf,0)) for nf in node_force]
                    
                else:
                    raise RuntimeError(f"Transition of type {trans_type} is unknown")
                    
                rate_data[f"Transition # {index} ({trans_type})"] = rate_trans
                    
        filename = os.path.join(pd.OUTPUT_DIR, f"rate_equations_iso_{i}.csv")    
        rate_data.to_csv(filename)
        
def calculate_rate_from_c_kinetics(c_kinetics, model_json_file):
    
    # Extract the rate laws from the kinetic data
       
    stretch = np.arange(X_MIN,X_MAX, X_STEP)    
    node_force = np.arange(F_MIN,F_MAX, F_STEP)
    
    with open(model_json_file, 'r') as f:
        mod = json.load(f)

    k_pc = mod["mybpc_parameters"]["c_k_stiff"]

    for i, iso in enumerate(c_kinetics):
        
        rate_data = pandas.DataFrame()
        
        rate_data["stretch"] = stretch
        
        rate_data["node_force"] = node_force
    
        for j, state in enumerate(iso):
            
            for new_state in state["transition"]:
                
                trans_type = new_state["rate_type"]
                trans_param = new_state["rate_parameters"]
                index = new_state["index"]
                
                if trans_type == "constant":
                    
                    rate_trans = [trans_param[0] for x in stretch]
                                   
                elif trans_type == "gaussian":
                    
                    rate_trans = [trans_param[0]*np.exp(-0.5 * k_pc * np.power(x + X_STEP/2, 2)/(1e18 * 1.38e-23*310)) for x in stretch]
                                    
                elif trans_type == "poly":
                    
                    rate_trans = [trans_param[0] + trans_param[1] * np.power(x, trans_param[2]) for x in stretch]
                    
                elif trans_type == "force_dependent":
                    
                    rate_trans = [trans_param[0] * (1 + trans_param[1] * max(nf,0)) for nf in node_force]
                    
                else:
                    raise RuntimeError(f"Transition of type {trans_type} is unknown")
                    
                rate_data[f"Transition # {index} ({trans_type})"] = rate_trans
                    
        filename = os.path.join(pd.OUTPUT_DIR, f"rate_equations_pc_iso_{i}.csv")    
        rate_data.to_csv(filename)