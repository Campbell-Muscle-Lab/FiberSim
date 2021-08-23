import os
import json

import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import pandas
import math

plt.rcParams.update({'font.family': "Arial"})
plt.rcParams.update({'font.size': 14}) 

# Add the half-sarcomere package

# ROOT = os.path.dirname(__file__)
# MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
# sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere
        
### Create stretch, node force and angle bins

# Stretch bins (nm)
    
X_MIN = -8
X_MAX = 8
X_STEP = 0.5

NB_INTER = int((X_MAX - X_MIN)/X_STEP)

# Force bin (nN)

F_MIN = -0.80
F_MAX = 0.80
F_STEP = 0.05 # same number of bin for stretch and node force

# Angle bins (°)
A_MIN = 90
A_MAX = 180
A_STEP = 9

NB_A_INTER = int((A_MAX - A_MIN)/A_STEP)

def compute_rate(model_file, protocol_file, dump_folder, output_folder, adj_bs = 0):
    """Approximates the governing rate functions for mybpc"""

    ### Get the hs_status dump files ###

    hs_file = []

    for filenames in os.listdir(dump_folder):
        filenames = os.path.join(dump_folder, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # sorting the dump files
    
    ### Get the time step ###

    protocol = np.loadtxt(protocol_file, skiprows=1)
    time_step = protocol[0][0]
            
    ### Get c_kinetics data

    # Extract the kinetics data
    
    c_kinetics = get_c_kinetics(model_file) 
    
    max_no_of_trans = c_kinetics[0][-1]["transition"][-1]["index"] + 1 
    
    attachement_trans = [] # List of indices for the 'a' type of transition
    
    ### Initialize transition arrays
    
    complete_transition = np.zeros((len(c_kinetics), max_no_of_trans, NB_INTER, NB_A_INTER),dtype=int)
    potential_transition = np.zeros((len(c_kinetics), max_no_of_trans, NB_INTER, NB_A_INTER),dtype=int)

    ### Fill the transition matrices
    
    # HS at t
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])
  
    for ind in range(1,len(hs_file)): 
        # HS at t + dt
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
      
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
                            pc_stretch_0 = get_stretch_interval(stretch)
                            
                            for j in range(0, 2*adj_bs+1): # Get the angle difference
                            
                                bs_ind = thick_fil_1[f"pc_nearest_a_n[x_{j}]"][pc_ind] # nearest BS have been recalculated
                                
                                found = False
                                if bs_ind == bs_ind_1: 
                                    found = True
                                    angle_diff = thick_fil_1[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                                    pc_align_factor_0 =  get_alignment_interval(angle_diff)
                                    
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
                                                
                                pc_stretch_0 = get_stretch_interval(stretch)
                                
                                angle_diff = thick_fil_1[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                                pc_align_factor_0 =  get_alignment_interval(angle_diff)
                                
                                complete_transition[pc_iso_0-1,idx,pc_stretch_0, pc_align_factor_0] += 1
                                                                
                    elif pc_0_type == 'A': # Get the pc stretch                    
                        
                        thin_ind = thick_fil_0["pc_bound_to_a_f"][pc_ind]
                        thin_fil_0 = hs_0["thin"][thin_ind]
                        bs_ind_0 = thick_fil_0["pc_bound_to_a_n"][pc_ind]
                        
                        # Find cb associated with this pc
                        pc_node_idx = thick_fil_0["pc_node_index"][pc_ind]  # crown index                           
                        cb_ind = pc_node_idx * 6    # crown index 
                    
                        stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil_0["bs_x"][bs_ind_0] 
                        pc_stretch_0 = get_stretch_interval(stretch)
                        
                        for j in range(0, 2*adj_bs+1): # Get the angle difference
                    
                            bs_ind = thick_fil_0[f"pc_nearest_a_n[x_{j}]"][pc_ind]  
                            found = False
                            
                            if bs_ind == bs_ind_0:
                                
                                found = True
                                angle_diff = thick_fil_0[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                                pc_align_factor_0 =  get_alignment_interval(angle_diff)
                                
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

                    pc_stretch_0 = get_stretch_interval(stretch)
                    
                    for j in range(0, 2*adj_bs+1): # Get the angle difference
            
                        bs_ind = thick_fil_0[f"pc_nearest_a_n[x_{j}]"][pc_ind]  
                        found = False
                        
                        if bs_ind == bs_ind_0:
                            found = True
                            angle_diff = thick_fil_0[f"pc_nearest_bs_angle_diff[x_{j}]"][pc_ind]                                    
                            pc_align_factor_0 =  get_alignment_interval(angle_diff)
                            
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
        
                                    pc_stretch_0 = get_stretch_interval(stretch[k])                                    
                                    pc_align_factor_0 =  get_alignment_interval(angle[k])
                                    potential_transition[pc_iso_0-1, idx_pot, pc_stretch_0, pc_align_factor_0] += 1 
                        
                        elif c_kinetics[pc_iso_0-1][new_state-1]["state_type"]== 'D': # Transition to another "D" state is always possible
                            for k in range(0, 2*adj_bs +1):   
                                pc_stretch_0 = get_stretch_interval(stretch[k])
                                pc_align_factor_0 =  get_alignment_interval(angle[k])
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
                
                    calculated_rate = np.zeros((NB_INTER, NB_A_INTER),dtype=float)
                    prob = np.zeros((NB_INTER, NB_A_INTER),dtype=float)
                    
                    p = np.zeros((NB_INTER, NB_A_INTER),dtype=float)
                    W = np.zeros((NB_INTER, NB_A_INTER),dtype=float)
                    
                    conf_interval_pos = np.zeros((NB_INTER, NB_A_INTER),dtype=float)
                    conf_interval_neg = np.zeros((NB_INTER, NB_A_INTER),dtype=float)
                                    
                    for stretch_bin in range(0,NB_INTER):
                        for angle_bin in range(0,NB_A_INTER):
                            
                                # Calculate rate
                            
                                prob[stretch_bin, angle_bin] = complete_transition[iso, trans, stretch_bin, angle_bin]/potential_transition[iso, trans, stretch_bin, angle_bin]
                                
                                if prob[stretch_bin, angle_bin] >= 1:
                                    print("Probability greater than 1")
                                    
                                calculated_rate[stretch_bin, angle_bin] = -np.log(1.0 - prob[stretch_bin, angle_bin]) / time_step # / -math.cos((angle[angle_bin] + A_STEP/2)*math.pi/180)

                                # Calculate CIs
                                
                                p[stretch_bin, angle_bin] = (complete_transition[iso, trans, stretch_bin, angle_bin] + 2)/(potential_transition[iso, trans, stretch_bin, angle_bin] + 4)
                                W[stretch_bin, angle_bin] = 2 * np.sqrt( (p[stretch_bin, angle_bin] * (1 - p[stretch_bin, angle_bin]))/ (potential_transition[iso, trans, stretch_bin, angle_bin] + 4) )
                                
                                conf_interval_pos[stretch_bin, angle_bin] = -np.log(1.0 - (p[stretch_bin, angle_bin] + W[stretch_bin, angle_bin]) ) / time_step #/ -math.cos((angle[angle_bin] + A_STEP/2)*math.pi/180)
                                conf_interval_neg[stretch_bin, angle_bin] = -np.log(1.0 - (p[stretch_bin, angle_bin] - W[stretch_bin, angle_bin]) ) / time_step #/ -math.cos((angle[angle_bin] + A_STEP/2)*math.pi/180)
    
                    trans_data["calculated_rate"] = calculated_rate.tolist()
                    trans_data["conf_pos"] = conf_interval_pos.tolist()
                    trans_data["conf_neg"] = conf_interval_neg.tolist()
                                    
                elif trans_type == 'd': # Detachment only depends on stretch
                
                    complete_transition_d = np.sum(complete_transition, axis = 3)
                    potential_transition_d = np.sum(potential_transition, axis = 3)
                
                    calculated_rate = np.zeros(NB_INTER,dtype=float)
                    prob = np.zeros(NB_INTER,dtype=float)
                    
                    p = np.zeros(NB_INTER,dtype=float)
                    W = np.zeros(NB_INTER,dtype=float)
                    
                    conf_interval_pos = np.zeros(NB_INTER, dtype=float)
                    conf_interval_neg = np.zeros(NB_INTER, dtype=float)
                
                    for stretch_bin in range(0,NB_INTER):
                        
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
        
    plot_rate(calculated_rates_dict, c_kinetics, model_file, output_folder)

def plot_rate(calculated_rates, c_kinetics, model_file, output_folder):
    
    angle = np.arange(A_MIN,A_MAX, A_STEP)

    for k, iso in enumerate(calculated_rates):
        
        rate_data = calculate_rate_from_c_kinetics(c_kinetics, model_file, k+1)
        
        for j, trans in enumerate(iso):
            
            output_image_file = os.path.join(output_folder, f"iso_{k}_transition_{j}.png")
            
            if trans["trans_type"] == 'a':
                
                idx = trans["trans_idx"]
                rate_type = trans["rate_type"]
                
                calc_rate = np.array(trans["calculated_rate"])
                conf_pos = np.array(trans["conf_pos"])
                conf_neg = np.array(trans["conf_neg"])
                
                # SUBPLOTS
                
                n_rows = int(np.floor(NB_A_INTER/2)) 
 
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
                                
                    align_factor = abs(-math.cos((angle[i] + A_STEP/2)*math.pi/180))
                    
                    true_rate = rate_data[f"Transition # {idx} ({rate_type})"] * align_factor
            
                    ax[i].errorbar(rate_data["stretch"], calc_rate[:,i], yerr=y_err, label = "calculated rate",
                                 ecolor = "tab:grey", fmt='-o', markersize = 2)

                    ax[i].plot(rate_data["stretch"], true_rate, label = f"rate law (angle diff = {angle[i] + A_STEP/2})")  
                    ax[0].set_title(f"Transition # {idx} ({rate_type})")
                    ax[i].legend(loc = "upper left")
                    ax[i].spines['right'].set_visible(False)
                    ax[i].spines['top'].set_visible(False)
                    ax[i].set_ylabel("Rate (s$^{-1}$)")
                    ax[i].set_ylim([0, max(rate_data[f"Transition # {idx} ({rate_type})"])])
                    
                ax[i].set_xlabel("CB stretch [nm]")
                
                # Second column
                
                for j in range(n_rows, n_rows + int(np.floor(NB_A_INTER/2))):  
                    
                    ax.append(fig.add_subplot(spec[j - n_rows, 1])) # add axes to second column
            
                    y_err = [
                      abs(calc_rate[:,j] - conf_neg[:,j]),
                      abs(calc_rate[:,j] - conf_pos[:,j])
                    ]
                            
                    align_factor = abs(-math.cos((angle[j] + A_STEP/2)*math.pi/180))
                    
                    ax[j].errorbar(rate_data["stretch"], calc_rate[:,j], yerr=y_err, label = "calculated rate",
                             ecolor = "tab:grey", fmt='-o', markersize = 2)
                    ax[j].set_ylabel("Rate (s$^{-1}$)")
                    ax[j].plot(rate_data["stretch"], rate_data[f"Transition # {idx} ({rate_type})"]*align_factor, label = f"rate law (angle diff = {angle[j] + A_STEP/2}")  
                    ax[n_rows].set_title(f"Transition # {idx} ({rate_type})")
                    ax[j].legend(loc = "upper left")
                    ax[j].spines['right'].set_visible(False)
                    ax[j].spines['top'].set_visible(False)
                    ax[j].set_ylim([0, max(rate_data[f"Transition # {idx} ({rate_type})"])])
                    
                ax[j].set_xlabel("CB stretch [nm]")
                plt.savefig(output_image_file, dpi = 150)
                    
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
                plt.errorbar(rate_data["stretch"] + X_STEP/2, calc_rate, yerr=y_err, label = "calculated rate",
                             ecolor = "tab:grey", fmt='-o', markersize = 2)
                plt.ylabel("Rate (s$^{-1}$)")
                plt.xlabel("CB stretch [nm]")
                plt.plot("stretch", f"Transition # {idx} ({rate_type})", data = rate_data, label = "rate law")  
                plt.title(f"Transition # {idx} ({rate_type})")
                plt.legend(loc = "upper left")
                plt.gca().spines['right'].set_visible(False)
                plt.gca().spines['top'].set_visible(False)
                plt.savefig(output_image_file, dpi = 150)
                                
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
                plt.ylabel("Rate (s$^{-1}$)")
                plt.xlim([0.5,1.5])
                plt.xticks([])
                
                true_rate = rate_data[f"Transition # {idx} ({rate_type})"].iloc[0]
                plt.ylim([true_rate - 1, true_rate + 1])
                plt.yticks([true_rate - 1, true_rate, true_rate +1])
                plt.plot([0,1,2], [true_rate, true_rate, true_rate], label = "rate law")   
                plt.title(f"Transition # {idx} ({rate_type})")
                plt.legend(loc = "upper left")
                plt.gca().spines['right'].set_visible(False)
                plt.gca().spines['top'].set_visible(False)
                plt.savefig(output_image_file, dpi = 150)
                
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
        
    return c_kinetics  
              
def calculate_rate_from_c_kinetics(c_kinetics, model_json_file, isotype = 1):
    
    # Extract the rate laws from the kinetic data
       
    stretch = np.arange(X_MIN,X_MAX, X_STEP)    
    node_force = np.arange(F_MIN,F_MAX, F_STEP)
    
    with open(model_json_file, 'r') as f:
        mod = json.load(f)

    k_pc = mod["mybpc_parameters"]["c_k_stiff"]
    
    rate_data = pandas.DataFrame()
    
    rate_data["stretch"] = stretch
    
    rate_data["node_force"] = node_force

    for j, state in enumerate(c_kinetics[isotype - 1]):
        
        for new_state in state["transition"]:
            
            trans_type = new_state["rate_type"]
            trans_param = new_state["rate_parameters"]
            index = new_state["index"]
            
            if trans_type == "constant":
                
                rate_trans = [trans_param[0] for x in stretch]
                               
            elif trans_type == "gaussian_pc":
                
                rate_trans = [trans_param[0]*np.exp(-0.5 * k_pc * np.power(x + X_STEP/2, 2)/(1e18 * 1.38e-23*310)) for x in stretch]
                                
            elif trans_type == "poly":
                
                rate_trans = [trans_param[0] + trans_param[1] * np.power(x, trans_param[2]) for x in stretch]
                
            elif trans_type == "force_dependent":
                
                rate_trans = [trans_param[0] * (1 + trans_param[1] * max(nf,0)) for nf in node_force]
                
            else:
                raise RuntimeError(f"Transition of type {trans_type} is unknown")
                
            rate_data[f"Transition # {index} ({trans_type})"] = rate_trans
                    
    return rate_data

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
