# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:08:32 2021

@author: sako231
"""

import os, sys
import json

ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

import path_definitions as pd
import kinetics_data as kd

import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt

import pandas
import math

def compute_m_kinetics_rate():
    """Approximates the governing rate functions of FiberSim"""

    ### Get the hs_status dump files ###

    hs_file = []

    for filenames in os.listdir(pd.HS_STATUS_FOLDER):
        filenames = os.path.join(pd.HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # sorting the dump files
    
    ### Get the time step ###

    protocol = np.loadtxt(pd.PROTOCOL_FILE, skiprows=1)
    time_step = protocol[0][0]
    
    ### Get adjacent_bs from option file
    
    with open(pd.OPTION_FILE, 'r') as f:
        opt = json.load(f)
        
    if 'adjacent_bs' in opt['options']:
        adj_bs = opt["options"]["adjacent_bs"]        
    else:
        adj_bs = 0
            
    ### Get m_kinetics data

    # Extract the kinetics data
    
    m_kinetics = kd.get_m_kinetics(pd.MODEL_FILE) 
    
    max_no_of_trans = m_kinetics[0][-1]["transition"][-1]["index"] + 1 
    
    ### Initialize transition arrays
    
    complete_transition = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=int)
    potential_transition = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=int)
    prob = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=float)
    calculated_rate = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=float)
               
    # HS at t
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])["hs_data"]
  
    for ind in range(1,len(hs_file)): 
        # HS at t + dt
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])["hs_data"]
      
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
                                                
                ### Determine if a transition occurred ###
                
                if cb_state_0 != cb_state_1: # A transition occurred
                
                    pot_trans = False
                   
                    # Find transition index
                    for trans in m_kinetics[cb_iso_0-1][cb_state_0-1]["transition"]:  # look through all transitions for cb_state_0
                        if trans["to"] == cb_state_1:
                            idx = trans["index"]
                            pot_trans = True
                            
                    if pot_trans == False:
                         raise RuntimeError(f"Transition index not found for transition from state {cb_state_0} to state {cb_state_1}")
                    
                    # Fill the element of the transition counter matrix
                    # with the proper transition index and stretch values
                    
                    if cb_0_type == 'S': 
                                         
                        if (cb_ind % 2) == 0: # Only even heads can "actively" transition
                                                 
                            # Get myosin node index
                            
                            node_index = int(np.floor(cb_ind/6))
                            
                            # Get node force
                            
                            node_force = thick_fil_0["node_forces"][node_index]
                            cb_node_force_0 = kd.get_node_force_interval(node_force)
                            
                            for j in range(0, 2*adj_bs+1): # Fill matrix for each angle
    
                                angle_diff = thick_fil_1[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    
                                cb_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                            
                                complete_transition[cb_iso_0-1,idx,cb_node_force_0, cb_align_factor_0] += 1
                            
                            # Check that dimer myosins transition together
                            if thick_fil_1["cb_state"][cb_ind] != thick_fil_1["cb_state"][cb_ind + 1]:
                                raise RuntimeError(f"Dimers did not follow the same transition from \
                                               state {cb_state_0} to state {cb_state_1}")                                                        
                                                       
                    elif cb_0_type == 'D':
                        
                        if cb_1_type == 'A': # Get the CB stretch
                            
                            thin_ind = thick_fil_1["cb_bound_to_a_f"][cb_ind]
                            thin_fil = hs_0["thin"][thin_ind] 
                            bs_ind_1 = thick_fil_1["cb_bound_to_a_n"][cb_ind]
                            
                            stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind_1] 
                            cb_stretch_0 = kd.get_stretch_interval(stretch)
                            
                            for j in range(0, 2*adj_bs+1): # Get the angle difference
                            
                                bs_ind = thick_fil_1[f"cb_nearest_a_n[x_{j}]"][cb_ind]    
                                
                                if bs_ind == bs_ind_1: 
                                    angle_diff = thick_fil_1[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    
                                    cb_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                                    complete_transition[cb_iso_0-1,idx,cb_stretch_0, cb_align_factor_0 ] += 1
                                    
                                # else:
                                #     print("not found")
                            
                        elif cb_1_type == 'D': 
                                                                      
                            thin_ind = thick_fil_0["cb_nearest_a_f"][cb_ind]
                            thin_fil = hs_0["thin"][thin_ind]
                            
                            for j in range(0, 2*adj_bs+1): # Loop over the nearest BS
            
                                bs_ind = thick_fil_1[f"cb_nearest_a_n[x_{j}]"][cb_ind]                                       
                                stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind] 
                                                
                                cb_stretch_0 = kd.get_stretch_interval(stretch)
                                
                                angle_diff = thick_fil_1[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    
                                cb_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                                
                                complete_transition[cb_iso_0-1,idx,cb_stretch_0, cb_align_factor_0] += 1
                                
                        elif cb_1_type == 'S':                        
                            
                            if (cb_ind % 2) == 0: # Only even heads can "actively" transition
                            
                                if thick_fil_0["cb_state"][cb_ind+1] == cb_state_0: # Transition can occur only if both dimer heads are in the same state at t
                        
                                    # Get myosin node index
                                
                                    node_index = int(np.floor(cb_ind/6))
                                
                                    # Get node force
                                    
                                    node_force = thick_fil_0["node_forces"][node_index]
                                    cb_node_force_0 = kd.get_node_force_interval(node_force)
                                    
                                    for j in range(0, 2*adj_bs+1): # Loop over angle difference
    
                                        angle_diff = thick_fil_1[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    
                                        cb_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                                    
                                        complete_transition[cb_iso_0-1,idx,cb_node_force_0,cb_align_factor_0] += 1
                                
                    elif cb_0_type == 'A': # Get the CB stretch                    
                        
                        thin_ind = thick_fil_0["cb_bound_to_a_f"][cb_ind]
                        thin_fil_0 = hs_0["thin"][thin_ind]
                        bs_ind_0 = thick_fil_0["cb_bound_to_a_n"][cb_ind]
                    
                        stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil_0["bs_x"][bs_ind_0] 
                        cb_stretch_0 = kd.get_stretch_interval(stretch)
                        
                        for j in range(0, 2*adj_bs+1): # Get the angle difference
                    
                            bs_ind = thick_fil_0[f"cb_nearest_a_n[x_{j}]"][cb_ind]  
                            
                            if bs_ind == bs_ind_0:
                                angle_diff = thick_fil_0[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    
                                cb_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                                complete_transition[cb_iso_0-1,idx,cb_stretch_0, cb_align_factor_0 ] += 1

                ### Determine all potential transitions depending on state type (S, D and A) ###
                
                # SRX (S)
                
                if cb_0_type == 'S':
                    
                    if (cb_ind % 2) == 0: # Only even heads can "actively" transition
                    
                        for trans in m_kinetics[cb_iso_0-1][cb_state_0-1]["transition"]: 

                            idx_pot = trans["index"]
                        
                            if thick_fil_0["cb_state"][cb_ind+1] == cb_state_0: # Transition can occur only if both dimer heads are in the same state at t
                         
                                # Get myosin node index
                                
                                node_index = int(np.floor(cb_ind/6))
                                
                                # Get node force
                                
                                node_force = thick_fil_0["node_forces"][node_index]
                                cb_node_force_0 = kd.get_node_force_interval(node_force)
                                
                                for j in range(0, 2*adj_bs+1): # Loop over angle differences
    
                                    angle = thick_fil_1[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    
                                    cb_align_factor_0 =  kd.get_alignment_interval(angle)
                               
                                    potential_transition[cb_iso_0-1,idx_pot,cb_node_force_0, cb_align_factor_0] += 1

                
                # Attached (A)
                
                elif cb_0_type == 'A': # Get CB stretch

                    thin_ind = thick_fil_0["cb_bound_to_a_f"][cb_ind]
                    bs_ind_0 = thick_fil_0["cb_bound_to_a_n"][cb_ind]
                    thin_fil = hs_0["thin"][thin_ind]   # Correction 
                    
                    stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind_0] 

                    cb_stretch_0 = kd.get_stretch_interval(stretch)
                    
                    for trans in m_kinetics[cb_iso_0-1][cb_state_0-1]["transition"]: 
                        
                        # All transitions from an attached state are possible
                        
                        idx_pot = trans["index"]
                        
                        found = False

                        for j in range(0, 2*adj_bs+1):
                    
                            bs_ind = thick_fil_0[f"cb_nearest_a_n[x_{j}]"][cb_ind]  
                            
                            if bs_ind == bs_ind_0:
                                angle = thick_fil_0[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    
                                cb_align_factor_0 =  kd.get_alignment_interval(angle)
                                potential_transition[cb_iso_0-1, idx_pot, cb_stretch_0, cb_align_factor_0] += 1 
                                found = True
                        if found == False:
                            print("not found")
                        
                # Detached (D)
                
                elif cb_0_type == 'D':
                    
                    # Get array of nearest BS and calculate stretch and angle difference array
                
                    bs_ind = np.zeros(2*adj_bs + 1, dtype = int)
                    stretch = np.zeros(2*adj_bs + 1)
                    angle = np.zeros(2*adj_bs + 1)
                    
                    thin_ind = thick_fil_0["cb_nearest_a_f"][cb_ind]
                    thin_fil = hs_0["thin"][thin_ind] # 
                    
                    for j in range(0, 2*adj_bs+1):
    
                        bs_ind[j] = thick_fil_1[f"cb_nearest_a_n[x_{j}]"][cb_ind] # nearest bs have been recalculated                                      
                        stretch[j] = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind[j]]
                        
                        angle[j] = thick_fil_1[f"cb_nearest_bs_angle_diff[x_{j}]"][cb_ind]                                    

                                           
                    for trans in m_kinetics[cb_iso_0-1][cb_state_0-1]["transition"]:
                        
                        new_state = trans["to"]
                        idx_pot = trans["index"]
                        
                        if m_kinetics[cb_iso_0-1][new_state-1]["state_type"]== 'A':
                            
                            # Check if the potential binding sites are available 
                            
                            for k in range(0, 2*adj_bs +1):
                                
                                bs_availability = thick_fil_1[f"cb_nearest_a_n_states[x_{k}]"][cb_ind] 
                            
                                if bs_availability == 2: # bs is free and available for attachment
        
                                    cb_stretch_0 = kd.get_stretch_interval(stretch[k])                                    
                                    cb_align_factor_0 =  kd.get_alignment_interval(angle[k])
                                    potential_transition[cb_iso_0-1, idx_pot, cb_stretch_0, cb_align_factor_0] += 1 

                        elif m_kinetics[cb_iso_0-1][new_state-1]["state_type"]== 'S':
                            
                                if (cb_ind % 2) == 0: # Only even heads can "actively" transition
                    
                                    if thick_fil_0["cb_state"][cb_ind+1] == cb_state_0: # Transition can occur only if both dimer heads are in the same state at t
                     
                                        # Get myosin node index
                                        
                                        node_index = int(np.floor(cb_ind/6))
                                        
                                        # Get node force
                                        
                                        node_force = thick_fil_0["node_forces"][node_index]
                                        cb_node_force_0 = kd.get_node_force_interval(node_force)
                                        
                                        for k in range(0, 2*adj_bs +1):
                                            cb_align_factor_0 =  kd.get_alignment_interval(angle[k])
                                            potential_transition[cb_iso_0-1,idx_pot,cb_node_force_0, cb_align_factor_0] += 1 
                        
                        elif m_kinetics[cb_iso_0-1][new_state-1]["state_type"]== 'D': # Transition to another "D" state is always possible
                            for k in range(0, 2*adj_bs +1):   
                                cb_stretch_0 = kd.get_stretch_interval(stretch[k])
                                cb_align_factor_0 =  kd.get_alignment_interval(angle[k])
                                potential_transition[cb_iso_0-1, idx_pot, cb_stretch_0, cb_align_factor_0[k]] += 1                                                           
                
                else:
                    raise RuntimeError(f"State #{cb_state_0} is neither type A, D or S")
                                
              
        hs_0 = hs_1 

    
    angle = [90, 108, 126, 144, 162]
    for iso in range(0, len(m_kinetics)):
        for trans in range(0, max_no_of_trans):
            for stretch_bin in range(0,kd.NB_INTER):
                for angle_bin in range(0,kd.NB_A_INTER):
                    if potential_transition[iso, trans, stretch_bin, angle_bin] == 0:
                        #print("No potential transition ! \n")
                        #print("# completed transition", complete_transition[iso, trans, stretch_bin])
                        prob[iso, trans, stretch_bin, angle_bin] = 0
                    else:
                        prob[iso, trans, stretch_bin, angle_bin] = complete_transition[iso, trans, stretch_bin, angle_bin]/potential_transition[iso, trans, stretch_bin, angle_bin]
                        if prob[iso, trans, stretch_bin, angle_bin] >= 1:
                            print("transition =", trans, "angle_bin =", angle_bin)
                            print("probability greater than 1")
                            prob[iso, trans, stretch_bin, angle_bin] = 1
    
                    calculated_rate[iso, trans, stretch_bin, angle_bin] = -np.log(1.0 - prob[iso, trans, stretch_bin, angle_bin]) / time_step / -math.cos((angle[angle_bin] + 18)*math.pi/180)
                           
    ### Calculate 95% CI ###
    
    p = (complete_transition + 2)/(potential_transition + 4)
    
    W = 2 * np.sqrt( (p * (1 - p))/ (potential_transition + 4) )

    conf_interval_pos = -np.log(1.0 - (p + W) ) / time_step    
    conf_interval_neg = -np.log(1.0 - (p - W) ) / time_step
   
    kd.calculate_rate_from_m_kinetics(m_kinetics, pd.MODEL_FILE) 
    
    for iso in range(0, len(m_kinetics)):

        rates_file = os.path.join(pd.OUTPUT_DIR, f"rate_equations_iso_{iso}.csv")  
        rate_data = pandas.read_csv(rates_file)
        
        for af in range(0, len(angle)):
        
            for j, col in enumerate(rate_data.columns):
                
                y_err = [
                  conf_interval_neg[iso, j-3, :],
                  conf_interval_pos[iso, j-3, :]
                ]
                
                filename = os.path.join(pd.OUTPUT_DIR, f"transition_{j-3}.png")
                            
                if "force_dependent" in col:
                    
                    plt.figure()
                    plt.plot(rate_data["node_force"], calculated_rate[iso, j-3, :, af], label = "calculated rate")
                    # plt.errorbar(rate_data["node_force"], calculated_rate[iso, j-3, :], yerr=y_err, label = "calculated rate", 
                    #              ecolor = "tab:grey", fmt='-o', errorevery = 2, markersize = 5)
                    plt.ylabel("Rate")
                    plt.xlabel("Node force (nN)")
                    plt.plot("node_force", col, data = rate_data, label = "rate law")  
                    plt.title(col)
                    plt.legend(loc = "upper left")
                    plt.savefig(filename)
                    plt.show()
                
                    
                elif "Transition" in col:
                    
                    plt.figure()
                    plt.plot(rate_data["stretch"], calculated_rate[iso, j-3, :, af], label = "calculated rate")
                    # plt.errorbar(rate_data["stretch"], calculated_rate[iso, j-3, :], yerr=y_err, label = "calculated rate", 
                    #              ecolor = "tab:grey", fmt='-o', errorevery = 2, markersize = 5)
                    plt.ylabel("Rate")
                    plt.xlabel("CB stretch [nm]")
                    plt.plot("stretch", col, data = rate_data, label = "rate law")  
                    plt.title(col)
                    plt.legend(loc = "upper left")
                    plt.savefig(filename)
                    plt.show()
                 
    
    return calculated_rate


calc_rate = compute_m_kinetics_rate()
