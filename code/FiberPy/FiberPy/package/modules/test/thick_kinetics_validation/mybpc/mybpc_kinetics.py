# -*- coding: utf-8 -*-
"""
Created on Thu May 13 22:13:24 2021

@author: srhko
"""

import os, sys
import json

ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", "..", "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

import path_definitions as pd
import kinetics_data as kd

import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

    ### Fill the transition matrices
    
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
                            
                            for trans in c_kinetics[pc_iso_0-1][pc_state_0-1]["transition"]: 
                                
                                # All transitions from an attached state are possible
                                
                                idx_pot = trans["index"]
                                
                                potential_transition[pc_iso_0-1, idx_pot, pc_stretch_0, pc_align_factor_0] += 1 
                                
                            break
                    
                    if found == False:
                        print("not found")
                        
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

    ### Calculate the rates
    
    calculated_rates_dict = []
    
    for iso in range(0, len(c_kinetics)): # Loop over the isotypes
    
        iso_data = []
        
        for i, state in enumerate(c_kinetics[iso]): # Loop over the states
            
            for j, transition in enumerate(state["transition"]): # For each state, loop over all possible transitions
                
                trans = transition["index"]
                trans_type = transition["trans_type"]
                rate_type = transition["rate_type"]
                
                trans_data = {}
                
                trans_data["trans_idx"] = trans
                trans_data["trans_type"] = trans_type
                trans_data["rate_type"] = rate_type

                if trans_type == 'a': # Attachment depends on stretch and angle
                
                    calculated_rate = np.zeros((kd.NB_INTER, kd.NB_A_INTER),dtype=float)
                    prob = np.zeros((kd.NB_INTER, kd.NB_A_INTER),dtype=float)
                    
                    p = np.zeros((kd.NB_INTER, kd.NB_A_INTER),dtype=float)
                    W = np.zeros((kd.NB_INTER, kd.NB_A_INTER),dtype=float)
                    
                    conf_interval_pos = np.zeros((kd.NB_INTER, kd.NB_A_INTER),dtype=float)
                    conf_interval_neg = np.zeros((kd.NB_INTER, kd.NB_A_INTER),dtype=float)
                                    
                    for stretch_bin in range(0,kd.NB_INTER):
                        for angle_bin in range(0,kd.NB_A_INTER):
                            
                                # Calculate rate
                            
                                prob[stretch_bin, angle_bin] = complete_transition[iso, trans, stretch_bin, angle_bin]/potential_transition[iso, trans, stretch_bin, angle_bin]
                                
                                if prob[stretch_bin, angle_bin] >= 1:
                                    print("Probability greater than 1")
                                    
                                calculated_rate[stretch_bin, angle_bin] = -np.log(1.0 - prob[stretch_bin, angle_bin]) / time_step # / -math.cos((angle[angle_bin] + kd.A_STEP/2)*math.pi/180)

                                # Calculate CIs
                                
                                p[stretch_bin, angle_bin] = (complete_transition[iso, trans, stretch_bin, angle_bin] + 2)/(potential_transition[iso, trans, stretch_bin, angle_bin] + 4)
                                W[stretch_bin, angle_bin] = 2 * np.sqrt( (p[stretch_bin, angle_bin] * (1 - p[stretch_bin, angle_bin]))/ (potential_transition[iso, trans, stretch_bin, angle_bin] + 4) )
                                
                                conf_interval_pos[stretch_bin, angle_bin] = -np.log(1.0 - (p[stretch_bin, angle_bin] + W[stretch_bin, angle_bin]) ) / time_step #/ -math.cos((angle[angle_bin] + kd.A_STEP/2)*math.pi/180)
                                conf_interval_neg[stretch_bin, angle_bin] = -np.log(1.0 - (p[stretch_bin, angle_bin] - W[stretch_bin, angle_bin]) ) / time_step #/ -math.cos((angle[angle_bin] + kd.A_STEP/2)*math.pi/180)
    
                    trans_data["calculated_rate"] = calculated_rate.tolist()
                    trans_data["conf_pos"] = conf_interval_pos.tolist()
                    trans_data["conf_neg"] = conf_interval_neg.tolist()
                                    
                elif trans_type == 'd': # Detachment only depends on stretch
                
                    complete_transition_d = np.sum(complete_transition, axis = 3)
                    potential_transition_d = np.sum(potential_transition, axis = 3)
                
                    calculated_rate = np.zeros(kd.NB_INTER,dtype=float)
                    prob = np.zeros(kd.NB_INTER,dtype=float)
                    
                    p = np.zeros(kd.NB_INTER,dtype=float)
                    W = np.zeros(kd.NB_INTER,dtype=float)
                    
                    conf_interval_pos = np.zeros(kd.NB_INTER, dtype=float)
                    conf_interval_neg = np.zeros(kd.NB_INTER, dtype=float)
                
                    for stretch_bin in range(0,kd.NB_INTER):
                        
                        # Calculate rate
                        
                        prob[stretch_bin] = complete_transition_d[iso, trans, stretch_bin]/potential_transition_d[iso, trans, stretch_bin]
                        
                        if prob[stretch_bin] >= 1:
                            print("Probability greater than 1")
                        calculated_rate[stretch_bin] = -np.log(1.0 - prob[stretch_bin]) / time_step
                        
                        # Calculate CIs
                                
                        p[stretch_bin] = (complete_transition_d[iso, trans, stretch_bin] + 2)/(potential_transition_d[iso, trans, stretch_bin] + 4)
                        W[stretch_bin] = 2 * np.sqrt( (p[stretch_bin] * (1 - p[stretch_bin]))/ (potential_transition_d[iso, trans, stretch_bin] + 4) )
                        
                        conf_interval_pos[stretch_bin] = -np.log(1.0 - (p[stretch_bin] + W[stretch_bin]) ) / time_step 
                        conf_interval_neg[stretch_bin] = -np.log(1.0 - (p[stretch_bin] - W[stretch_bin]) ) / time_step 

                    trans_data["calculated_rate"] = calculated_rate.tolist()
                    trans_data["conf_pos"] = conf_interval_pos.tolist()
                    trans_data["conf_neg"] = conf_interval_neg.tolist()                    

                elif trans_type == 'x': # Transition with a constant rate law
                
                    complete_transition_x = np.sum(np.sum(complete_transition, axis = 3), axis = 2)
                    potential_transition_x = np.sum(np.sum(potential_transition, axis = 3), axis = 2)
                
                    calculated_rate = 0.0
                    prob = 0.0
                    
                    p = 0.0
                    W = 0.0
                    
                    conf_interval_pos = 0.0
                    conf_interval_neg = 0.0
                    
                    # Calculate rate

                    prob = complete_transition_x[iso, trans]/potential_transition_x[iso, trans]
                    
                    if prob >= 1:
                        print("Probability greater than 1")
                        
                    calculated_rate = -np.log(1.0 - prob) / time_step
                    
                    # Calculate CIs
                            
                    p = (complete_transition_x[iso, trans] + 2)/(potential_transition_x[iso, trans] + 4)
                    W = 2 * np.sqrt( (p * (1 - p))/ (potential_transition_x[iso, trans] + 4) )
                    
                    conf_interval_pos = -np.log(1.0 - (p + W) ) / time_step 
                    conf_interval_neg = -np.log(1.0 - (p - W) ) / time_step 
    
                    trans_data["calculated_rate"] = calculated_rate.tolist()
                    trans_data["conf_pos"] = conf_interval_pos.tolist()
                    trans_data["conf_neg"] = conf_interval_neg.tolist()
                    
                iso_data.append(trans_data)
                
        calculated_rates_dict.append(iso_data)
                
    # Save the calculated rates structure in a JSON file
    outfile = os.path.join(pd.OUTPUT_DIR, "calc_dict.json")  
    with open(outfile, 'w') as f:
        json.dump(calculated_rates_dict, f)

    # Load calculated rates

    outfile = os.path.join(pd.OUTPUT_DIR, "calc_dict.json")  
    with open(outfile, 'r') as f:
        calculated_rates_dict = json.load(f)
        
    ### PLOTS
    
    angle = np.arange(kd.A_MIN,kd.A_MAX, kd.A_STEP)

    for i, iso in enumerate(calculated_rates_dict):
        
        rates_file = os.path.join(pd.OUTPUT_DIR, f"rate_equations_iso_{i}.csv")  
        rate_data = pandas.read_csv(rates_file)
        
        for j, trans in enumerate(iso):
            
            if trans["trans_type"] == 'a':
                
                idx = trans["trans_idx"]
                rate_type = trans["rate_type"]
                
                calc_rate = np.array(trans["calculated_rate"])
                conf_pos = np.array(trans["conf_pos"])
                conf_neg = np.array(trans["conf_neg"])
                
                # SUBPLOTS
                
                n_rows = int(np.floor(kd.NB_A_INTER/2)) 
 
                # Now create figure
                fig = plt.figure(constrained_layout=False, figsize = [15,15])
                spec = gridspec.GridSpec(nrows=n_rows,
                                         ncols=2,
                                         figure=fig)
                ax = []
                
                # First column
                
                for i in range(0, n_rows):
                    
                    ax.append(fig.add_subplot(spec[i, 0])) # add axes to first column
            
                    y_err = [
                      abs(calc_rate[:,i] - conf_neg[:,i]),
                      abs(calc_rate[:,i] - conf_pos[:,i])
                    ]
                                
                    align_factor = abs(-math.cos((angle[i] + kd.A_STEP/2)*math.pi/180))
                    
                    true_rate = rate_data[f"Transition # {idx} ({rate_type})"] * align_factor
            
                    ax[i].errorbar(rate_data["stretch"], calc_rate[:,i], yerr=y_err, label = "calculated rate",
                                 ecolor = "tab:grey", fmt='-o', markersize = 2)

                    ax[i].plot(rate_data["stretch"], true_rate, label = f"rate law (angle diff = {angle[i] + kd.A_STEP/2})")  
                    ax[0].set_title(f"Transition # {idx} ({rate_type})")
                    ax[i].legend(loc = "upper left")
                    ax[i].spines['right'].set_visible(False)
                    ax[i].spines['top'].set_visible(False)
                    ax[i].set_ylabel("Rate")
                    ax[i].set_ylim([0, max(rate_data[f"Transition # {idx} ({rate_type})"])])
                    
                ax[i].set_xlabel("CB stretch [nm]")
                
                # Second column
                
                for j in range(n_rows, n_rows + int(np.floor(kd.NB_A_INTER/2))):  
                    
                    ax.append(fig.add_subplot(spec[j - n_rows, 1])) # add axes to second column
            
                    y_err = [
                      abs(calc_rate[:,j] - conf_neg[:,j]),
                      abs(calc_rate[:,j] - conf_pos[:,j])
                    ]
                            
                    align_factor = abs(-math.cos((angle[j] + kd.A_STEP/2)*math.pi/180))
                    
                    ax[j].errorbar(rate_data["stretch"], calc_rate[:,j], yerr=y_err, label = "calculated rate",
                             ecolor = "tab:grey", fmt='-o', markersize = 2)
                    ax[j].set_ylabel("Rate")
                    ax[j].plot(rate_data["stretch"], rate_data[f"Transition # {idx} ({rate_type})"]*align_factor, label = f"rate law (angle diff = {angle[j] + kd.A_STEP/2}")  
                    ax[n_rows].set_title(f"Transition # {idx} ({rate_type})")
                    ax[j].legend(loc = "upper left")
                    ax[j].spines['right'].set_visible(False)
                    ax[j].spines['top'].set_visible(False)
                    ax[j].set_ylim([0, max(rate_data[f"Transition # {idx} ({rate_type})"])])
                    
                ax[j].set_xlabel("CB stretch [nm]")
                plt.show()
                    
            elif trans["trans_type"] == 'd':
                
                idx = trans["trans_idx"]
                rate_type = trans["rate_type"]
                calc_rate = np.array(trans["calculated_rate"])
                conf_pos = np.array(trans["conf_pos"])
                conf_neg = np.array(trans["conf_neg"])
                
                y_err = [
                      abs(calc_rate - conf_neg),
                      abs(calc_rate - conf_pos)
                    ]
                
                plt.figure()
                plt.errorbar(rate_data["stretch"] + kd.X_STEP/2, calc_rate, yerr=y_err, label = "calculated rate",
                             ecolor = "tab:grey", fmt='-o', markersize = 2)
                plt.ylabel("Rate")
                plt.xlabel("CB stretch [nm]")
                plt.plot("stretch", f"Transition # {idx} ({rate_type})", data = rate_data, label = "rate law")  
                plt.title(f"Transition # {idx} ({rate_type})")
                plt.legend(loc = "upper left")
                plt.show()
                                
            elif trans["trans_type"] == 'x':
                
                idx = trans["trans_idx"]
                rate_type = trans["rate_type"]
                calc_rate = np.array(trans["calculated_rate"])
                conf_pos = np.array(trans["conf_pos"])
                conf_neg = np.array(trans["conf_neg"])
                
                y_err = [
                      abs(calc_rate - conf_neg),
                      abs(calc_rate - conf_pos)
                    ]
                
                plt.figure()
                plt.errorbar([1,2], [calc_rate, calc_rate], yerr=y_err, label = "calculated rate",
                             ecolor = "tab:grey", fmt='o', markersize = 2)
                plt.ylabel("Rate")
                plt.xlim([0.5,1.5])
                plt.xticks([])
                
                true_rate = rate_data[f"Transition # {idx} ({rate_type})"].iloc[0]
                plt.ylim([true_rate - 1, true_rate + 1])
                plt.yticks([true_rate - 1, true_rate, true_rate +1])
                plt.plot([0,1,2], [true_rate, true_rate, true_rate], label = "rate law")   
                plt.title(f"Transition # {idx} ({rate_type})")
                plt.legend(loc = "upper left")
                plt.show()

compute_c_kinetics_rate()
