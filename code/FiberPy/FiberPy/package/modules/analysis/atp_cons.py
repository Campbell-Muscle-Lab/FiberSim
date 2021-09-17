# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 14:11:32 2021

@author: sako231
"""

import os
import json

import numpy as np
import pandas as pd
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# Add the half-sarcomere package

# ROOT = os.path.dirname(__file__)
# MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
# sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

def get_ATP_cons(data, batch_file_string): # Calculate the ATP consumption rate 
    
    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)
    
    model_array = []
    dump_array = []
    
    for elmt in data["model_file"]:

        if (data['relative_to'] == 'this_file'):
            model_file = os.path.join(base_folder, elmt)
        else:
            model_file = elmt
    
        model_array.append(model_file)
    
    for elmt in data["status_folder"]:

        if (data['relative_to'] == 'this_file'):
            dump_folder = os.path.join(base_folder, elmt)
        else:
            dump_folder = elmt    
        dump_array.append(dump_folder)
    
    # Create lists to hold data
    case = []
    ATP_cons_per_head = []
    ATP_cons = []
    
    for i in range(len(model_array)):
        
        ATP_head, ATP_cons_mol = calculate_ATPase_rate(model_array[i], dump_array[i])
        
        case.append(i+1)
        ATP_cons_per_head.append(ATP_head)
        ATP_cons.append(ATP_cons_mol)
    
    # If more than one rate is calculated, create a figure
    ### CREATE figure

    # Take list and create data frame
    
    pd.set_option('max_colwidth', 400)

    r = pd.DataFrame({'case': case,
                    'ATP consumption \n rate \n (# molecule / myosin head /s)': ATP_cons_per_head,
                    'ATP \n consumption \n rate \n (mol / s / mg)': ATP_cons})

    # Save the data as an excel file if required in the batch file
    if(data['output_data_file_string']):
        if (data['relative_to'] == 'this_file'):
            output_file_string = os.path.join(base_folder,
                                              data['output_data_file_string'])
        else:
            output_file_string = data['output_data_file_string']

        print('Writing ATP consumption rates to %s' % output_file_string)
        
        
        #writer = pd.ExcelWriter(output_file_string) 
        
        r.to_excel(output_file_string, sheet_name='atp_consumption',
                   engine='openpyxl',
                   index=False)
        

               

def calculate_ATPase_rate(model_file, dump_folder):
    
    """Calculate the number of ATP molecules consumed per unit time"""

    ### Get the hs_status dump files ###

    hs_file = []

    for filenames in os.listdir(dump_folder):
        filenames = os.path.join(dump_folder, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # sorting the dump files    
    
    ### Initialize the number of ATPs consumed
    
    no_of_ATP = 0
    
    ### Extract the kinetics data
    
    m_kinetics = get_m_kinetics(model_file) 
    
    ### Loop over the dump files
               
    # HS at t
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])
    
    # Get the extension array
    
    ext = hs_0["hs_data"]["cb_extensions"]
  
    for ind in range(1,len(hs_file)): 
        # HS at t + dt
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
      
        for i, thick_fil_0 in enumerate(hs_0["thick"]): # Loop over all thick filaments
            
            thick_fil_1 = hs_1["thick"][i] 
   
            for cb_ind, state_0 in enumerate(thick_fil_0["cb_state"]): # loop over all CBs
                
                cb_iso_0 = thick_fil_0["cb_iso"][cb_ind] # CB isotype at t
                cb_iso_1 = thick_fil_1["cb_iso"][cb_ind] # CB isotype at t + dt
               
                # Check isotype does not change through time
                if cb_iso_0 != cb_iso_1:                    
                    raise RuntimeError(f"Isotype #{cb_iso_0} turned into isotype #{cb_iso_1}")
                    
                cb_state_0 = thick_fil_0["cb_state"][cb_ind] # CB state at t
                cb_state_1 = thick_fil_1["cb_state"][cb_ind] # CB state at t + dt
                
                cb_0_type = m_kinetics[cb_iso_0-1][cb_state_0-1]["state_type"] # CB type at t
                cb_1_type = m_kinetics[cb_iso_1-1][cb_state_1-1]["state_type"] # CB type at t + dt
                                                
                ### Determine if a detachment occurred ### 
           
                if cb_0_type == 'A' and cb_1_type == 'D': # A detachment event occured 
                    
                    cb_ext = ext[cb_state_0  - 1] # get the CB extension
                    
                    if cb_ext > 0.0 : # A post-power stroke detachment event occured 
                    
                        no_of_ATP += 1                                
              
        hs_0 = hs_1
        
    total_time = hs_0["hs_data"]["time"] # Time in s

    no_myosin_heads = len(hs_0["thick"]) * hs_0["hs_data"]["m_nodes_per_thick_filament"] * hs_0 ["thick"][0]["m_cbs_per_node"]
    
    ATP_per_head = no_of_ATP/no_myosin_heads/total_time # Get the number of ATP molecules consumed per unit time per myosin head
    
    N_A = 6e23 # Avogadro number
    
    vol = len(hs_0["thick"])* hs_0["hs_data"]["hs_length"] *1e-9 /hs_0["hs_data"]["m_filament_density"] # half-sarcomere volume (in m³)
    
    vol = vol * 1e3 # volume in liter

    density = 1.055 # density in g/ml (PMC6461811)
    

    ATP_cons = no_of_ATP / N_A # get mol
    
    ATP_cons = ATP_cons / total_time # get mol/s

    ATP_cons = ATP_cons / (density * vol) # get mol/s/kg

    ATP_cons = ATP_cons * 1e-6 # get mol/s/mg
                
    return ATP_per_head, ATP_cons
    
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
                
    return m_kinetics
