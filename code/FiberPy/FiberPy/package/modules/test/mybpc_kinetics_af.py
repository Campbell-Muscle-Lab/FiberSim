# -*- coding: utf-8 -*-
"""
Created on Thu May 13 22:13:24 2021

@author: srhko
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

def compute_c_kinetics_rate():
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
            
    ### Get c_kinetics data

    # Extract the kinetics data
    
    c_kinetics = kd.get_c_kinetics(pd.MODEL_FILE) 
    
    max_no_of_trans = c_kinetics[0][-1]["transition"][-1]["index"] + 1 
    
    attachement_trans = [] # List of indices for the 'a' type of transition
    
    ### Initialize transition arrays
    
    complete_transition = np.zeros((len(c_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=int)
    potential_transition = np.zeros((len(c_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=int)
    prob = np.zeros((len(c_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=float)
    calculated_rate = np.zeros((len(c_kinetics), max_no_of_trans, kd.NB_INTER, kd.NB_A_INTER),dtype=float)
               
    # HS at t
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])["hs_data"]
  
    for ind in range(1,len(hs_file)): 
        # HS at t + dt
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])["hs_data"]
      
        for i, thick_fil_0 in enumerate(hs_0["thick"]): # Loop over all thick filaments
            
            thick_fil_1 = hs_1["thick"][i] 
   
            for pc_ind, state_0 in enumerate(thick_fil_0["pc_state"]): # loop over all pcs
                
                pc_iso_0 = thick_fil_0["pc_iso"][pc_ind] # pc isotype at t
                pc_iso_1 = thick_fil_1["pc_iso"][pc_ind] # pc isotype at t + dt
               
                # Check isotype does not change through time
                if pc_iso_0 != pc_iso_1:                    
                    raise RuntimeError(f"Isotype #{pc_iso_0} turned into isotype #{pc_iso_1}")
                    
                pc_state_0 = thick_fil_0["pc_state"][pc_ind] # CB state at t
                pc_state_1 = thick_fil_1["pc_state"][pc_ind] # CB state at t + dt
                
                pc_0_type = c_kinetics[pc_iso_0-1][pc_state_0-1]["state_type"] # CB type at t
                pc_1_type = c_kinetics[pc_iso_1-1][pc_state_1-1]["state_type"] # CB type at t + dt
                                                
                ### Determine if a transition occurred ###
                
                if pc_state_0 != pc_state_1: # A transition occurred
                
                    pot_trans = False
                   
                    # Find transition index
                    for trans in c_kinetics[pc_iso_0-1][pc_state_0-1]["transition"]:  # look through all transitions for pc_state_0
                        if trans["to"] == pc_state_1:
                            idx = trans["index"]
                            pot_trans = True
                            
                    if pot_trans == False:
                         raise RuntimeError(f"Transition index not found for transition from state {pc_state_0} to state {pc_state_1}")
                    
                    # Fill the element of the transition counter matrix
                    # with the proper transition index and stretch values
                    
                    if pc_0_type == 'D':
                        
                        if pc_1_type == 'A': # Get the pc stretch
                        
                            if idx not in attachement_trans: 
                                attachement_trans.append(idx) 
                            
                            thin_ind = thick_fil_1["pc_bound_to_a_f"][pc_ind]
                            thin_fil = hs_0["thin"][thin_ind] 
                            bs_ind_1 = thick_fil_1["pc_bound_to_a_n"][pc_ind]
                            
                            # Find cb associated with this pc
                            pc_node_idx = thick_fil_0["pc_node_index"][pc_ind] # crown index                           
                            cb_ind = pc_node_idx * 6 # cb index
                            
                            stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind_1] # x at t0 !
                            pc_stretch_0 = kd.get_stretch_interval(stretch)
                            
                            for j in range(0, 2*adj_bs+1): # Get the angle difference
                            
                                bs_ind = thick_fil_1[f"pc_nearest_a_n[x_{j}]"][pc_ind] # nearest BS have been recalculated
                                
                                found = False
                                if bs_ind == bs_ind_1: 
                                    found = True
                                    angle_diff = thick_fil_1[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                                    pc_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                                    
                                    complete_transition[pc_iso_0-1,idx,pc_stretch_0, pc_align_factor_0 ] += 1
                                    break
                                
                            if found == False:
                                print("not found")
                            
                        elif pc_1_type == 'D': 
                                                                      
                            thin_ind = thick_fil_0["pc_nearest_a_f"][pc_ind]
                            thin_fil = hs_0["thin"][thin_ind]
                            
                            for j in range(0, 2*adj_bs+1): # Loop over the nearest BS
            
                                bs_ind = thick_fil_1[f"pc_nearest_a_n[x_{j}]"][pc_ind] 

                                # Find cb associated with this pc
                                pc_node_idx = thick_fil_0["pc_node_index"][pc_ind]  # crown index                          
                                cb_ind = pc_node_idx * 6
                                      
                                stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind] # x at t0 !
                                                
                                pc_stretch_0 = kd.get_stretch_interval(stretch)
                                
                                angle_diff = thick_fil_1[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                                pc_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                                
                                complete_transition[pc_iso_0-1,idx,pc_stretch_0, pc_align_factor_0] += 1
                                                                
                    elif pc_0_type == 'A': # Get the pc stretch                    
                        
                        thin_ind = thick_fil_0["pc_bound_to_a_f"][pc_ind]
                        thin_fil_0 = hs_0["thin"][thin_ind]
                        bs_ind_0 = thick_fil_0["pc_bound_to_a_n"][pc_ind]
                        
                        # Find cb associated with this pc
                        pc_node_idx = thick_fil_0["pc_node_index"][pc_ind]  # crown index                           
                        cb_ind = pc_node_idx * 6    # crown index 
                    
                        stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil_0["bs_x"][bs_ind_0] 
                        pc_stretch_0 = kd.get_stretch_interval(stretch)
                        
                        for j in range(0, 2*adj_bs+1): # Get the angle difference
                    
                            bs_ind = thick_fil_0[f"pc_nearest_a_n[x_{j}]"][pc_ind]  
                            found = False
                            
                            if bs_ind == bs_ind_0:
                                
                                found = True
                                angle_diff = thick_fil_0[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                                pc_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                                
                                complete_transition[pc_iso_0-1,idx,pc_stretch_0, pc_align_factor_0 ] += 1
                                break
                            
                        if found == False:
                            print("not found")
                                
                ### Determine all potential transitions depending on state type (D and A) ###
                                
                # Attached (A)
                
                if pc_0_type == 'A': # Get pc stretch

                    thin_ind = thick_fil_0["pc_bound_to_a_f"][pc_ind]
                    bs_ind_0 = thick_fil_0["pc_bound_to_a_n"][pc_ind]
                    thin_fil = hs_0["thin"][thin_ind]  
                    
                    # Find cb associated with this pc
                    pc_node_idx = thick_fil_1["pc_node_index"][pc_ind]  # crown index                           
                    cb_ind = pc_node_idx * 6
                    
                    stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind_0] 

                    pc_stretch_0 = kd.get_stretch_interval(stretch)
                    
                    for j in range(0, 2*adj_bs+1): # Get the angle difference
            
                        bs_ind = thick_fil_0[f"pc_nearest_a_n[x_{j}]"][pc_ind]  
                        found = False
                        
                        if bs_ind == bs_ind_0:
                            found = True
                            angle_diff = thick_fil_0[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                            pc_align_factor_0 =  kd.get_alignment_interval(angle_diff)
                            break
                    
                    if found == False:
                        print("not found")
                    
                    for trans in c_kinetics[pc_iso_0-1][pc_state_0-1]["transition"]: 
                        
                        # All transitions from an attached state are possible
                        
                        idx_pot = trans["index"]
                        
                        potential_transition[pc_iso_0-1, idx_pot, pc_stretch_0, pc_align_factor_0] += 1 
                        
                # Detached (D)
                
                elif pc_0_type == 'D':
                    
                    # Get array of nearest BS and calculate stretch and angle difference array
                
                    bs_ind = np.zeros(2*adj_bs + 1, dtype = int)
                    stretch = np.zeros(2*adj_bs + 1)
                    angle = np.zeros(2*adj_bs + 1)
                    
                    thin_ind = thick_fil_0["pc_nearest_a_f"][pc_ind]
                    thin_fil = hs_0["thin"][thin_ind] 
                    
                    for j in range(0, 2*adj_bs+1):
    
                        bs_ind[j] = thick_fil_1[f"pc_nearest_a_n[x_{j}]"][pc_ind] # nearest bs have been recalculated                                      
                        
                        # Find cb associated with this pc
                        pc_node_idx = thick_fil_1["pc_node_index"][pc_ind] # crown index                            
                        cb_ind = pc_node_idx * 6
                        
                        stretch[j] = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind[j]] # x has not been recalculated
                        
                        angle[j] = thick_fil_1[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]  # nearest bs have been recalculated                                    

                                           
                    for trans in c_kinetics[pc_iso_0-1][pc_state_0-1]["transition"]:
                        
                        new_state = trans["to"]
                        idx_pot = trans["index"]
                        
                        if c_kinetics[pc_iso_0-1][new_state-1]["state_type"]== 'A':
                            
                            # Check if the potential binding sites are available 
                            
                            for k in range(0, 2*adj_bs +1):
                                
                                bs_availability = thick_fil_1[f"pc_nearest_a_n_states[x_{k}]"][pc_ind] 
                            
                                if bs_availability == 2: # bs is free and available for attachment
        
                                    pc_stretch_0 = kd.get_stretch_interval(stretch[k])                                    
                                    pc_align_factor_0 =  kd.get_alignment_interval(angle[k])
                                    potential_transition[pc_iso_0-1, idx_pot, pc_stretch_0, pc_align_factor_0] += 1 
                        
                        elif c_kinetics[pc_iso_0-1][new_state-1]["state_type"]== 'D': # Transition to another "D" state is always possible
                            for k in range(0, 2*adj_bs +1):   
                                pc_stretch_0 = kd.get_stretch_interval(stretch[k])
                                pc_align_factor_0 =  kd.get_alignment_interval(angle[k])
                                potential_transition[pc_iso_0-1, idx_pot, pc_stretch_0, pc_align_factor_0[k]] += 1                                                           
                
                else:
                    raise RuntimeError(f"State #{pc_state_0} is neither type A, D or S")
                                
              
        hs_0 = hs_1 

    angle = np.arange(kd.A_MIN,kd.A_MAX, kd.A_STEP)
    for iso in range(0, len(c_kinetics)):
        for trans in range(0, max_no_of_trans):
            for stretch_bin in range(0,kd.NB_INTER):
                for angle_bin in range(0,kd.NB_A_INTER):
                    if potential_transition[iso, trans, stretch_bin, angle_bin] == 0:
                        #print("No potential transition ! \n")
                        prob[iso, trans, stretch_bin, angle_bin] = 0
                    else:
                        prob[iso, trans, stretch_bin, angle_bin] = complete_transition[iso, trans, stretch_bin, angle_bin]/potential_transition[iso, trans, stretch_bin, angle_bin]
                        if prob[iso, trans, stretch_bin, angle_bin] >= 1:
                            print("probability greater than 1")
                            prob[iso, trans, stretch_bin, angle_bin] = 1
                    if trans in attachement_trans:
                        calculated_rate[iso, trans, stretch_bin, angle_bin] = -np.log(1.0 - prob[iso, trans, stretch_bin, angle_bin]) / time_step / -math.cos((angle[angle_bin] + 18)*math.pi/180)
                    else:
                        calculated_rate[iso, trans, stretch_bin, angle_bin] = -np.log(1.0 - prob[iso, trans, stretch_bin, angle_bin]) / time_step

                        
    ### Calculate 95% CI ###
    
    p = (complete_transition + 2)/(potential_transition + 4)
    
    W = 2 * np.sqrt( (p * (1 - p))/ (potential_transition + 4) )

    conf_interval_pos = -np.log(1.0 - (p + W) ) / time_step    
    conf_interval_neg = -np.log(1.0 - (p - W) ) / time_step
   
    kd.calculate_rate_from_c_kinetics(c_kinetics, pd.MODEL_FILE) 
    
    for iso in range(0, len(c_kinetics)):

        rates_file = os.path.join(pd.OUTPUT_DIR, f"rate_equations_pc_iso_{iso}.csv")  
        rate_data = pandas.read_csv(rates_file)
        
        for af in range(0, len(angle)):
        
            for j, col in enumerate(rate_data.columns):
                
                # y_err = [
                #   conf_interval_neg[iso, j-3, :],
                #   conf_interval_pos[iso, j-3, :]
                # ]
                
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


calc_rate = compute_c_kinetics_rate()
