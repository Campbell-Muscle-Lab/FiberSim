# -*- coding: utf-8 -*-
"""
Created on Mon May  3 17:48:49 2021

@author: sako231
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 20:51:24 2021

@author: srhko
"""

import os
import json

from modules.half_sarcomere import half_sarcomere

import path_definitions as pd
import kinetics_data as kd

import numpy as np
from natsort import os_sorted
import matplotlib.pyplot as plt

def compute_m_kinetics_rate():
    """Approximates the governing rate functions of FiberSim"""

    ### Get the hs_status dump files ###

    hs_file = []

    for filenames in os.listdir(pd.HS_STATUS_FOLDER):
        filenames = os.path.join(pd.HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs_file = os_sorted(hs_file) # sorting the dump files
    
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
    
    complete_transition = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER),dtype=int)
    potential_transition = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER),dtype=int)
    prob = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER),dtype=float)
    calculated_rate = np.zeros((len(m_kinetics), max_no_of_trans, kd.NB_INTER),dtype=float)
               
    # HS at t0
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])["hs_data"]
  
    for ind in range(1,len(hs_file)):
                
        # HS at t0 + delta t
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])["hs_data"]
      
        for i, thick_fil_0 in enumerate(hs_0["thick"]):
            
            thick_fil_1 = hs_1["thick"][i]
   
            for cb_ind, state_0 in enumerate(thick_fil_0["cb_state"]):
                
                cb_iso_0 = thick_fil_0["cb_iso"][cb_ind]
                cb_iso_1 = thick_fil_1["cb_iso"][cb_ind]
               
                # Check isotype does not change through time
                if cb_iso_0 != cb_iso_1:                    
                    raise RuntimeError(f"Isotype #{cb_iso_0} turned into isotype #{cb_iso_1}")
                    
                cb_state_0 = thick_fil_0["cb_state"][cb_ind] # CB state at t0
                cb_state_1 = thick_fil_1["cb_state"][cb_ind] # CB state at t0 + time_step
                
                cb_0_type = m_kinetics[cb_iso_0-1][cb_state_0-1]["state_type"] # CB type at t0 + time_step
                cb_1_type = m_kinetics[cb_iso_0-1][cb_state_1-1]["state_type"] # CB type at t0 + time_step
                                                
                ### Determine if a transition occurred ###
                
                if cb_state_0 != cb_state_1: # A transition occurred
                    
                    # Find transition index
                    for trans in m_kinetics[cb_iso_0-1][cb_state_0-1]["transition"]:  # look through all transitions for cb_state_0
                        if trans["to"] == cb_state_1:
                            idx = trans["index"]
                    
                    # Fill the element of the transition counter matrix
                    # with the proper transition index and stretch values
                    
                    if cb_0_type == 'S':
                    
                        # Check is cb_ind is odd -> no complete transition to account for
                        # Complete transition if cb_ind is even
                        # Check if the SRX dimers transition at the same time
                        break
                                      
                    elif cb_0_type == 'A': # Get the CB stretch
                        
                        thin_ind = thick_fil_0["cb_bound_to_a_f"][cb_ind]
                        thin_fil = hs_0["thin"][thin_ind]
                        bs_ind_0 = thick_fil_0["cb_bound_to_a_n"][cb_ind]
                    
                        stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind_0]                     
                        cb_stretch_0 = kd.get_stretch_interval(stretch)
                                            
                        complete_transition[cb_iso_0-1,idx,cb_stretch_0] += 1
                        
                    elif cb_0_type == 'D':
                        
                        if cb_1_type == 'A': # Get the CB stretch
                            
                            thin_ind = thick_fil_1["cb_bound_to_a_f"][cb_ind]
                            thin_fil = hs_0["thin"][thin_ind]
                            bs_ind_1 = thick_fil_1["cb_bound_to_a_n"][cb_ind]
                            
                            stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind_1]                     
                            cb_stretch_0 = kd.get_stretch_interval(stretch)
                                            
                            complete_transition[cb_iso_0-1,idx,cb_stretch_0] += 1
                            
                        elif cb_1_type == 'D': 
                            
                        # Get array of nearest BS and calculate stretch array
                
                            bs_ind = np.zeros(2*adj_bs + 1, dtype = int)
                            stretch = np.zeros(2*adj_bs + 1)
                            
                            thin_ind = thick_fil_0["cb_nearest_a_f"][cb_ind]
                            thin_fil = hs_0["thin"][thin_ind]
                            
                            for j in range(0, 2*adj_bs+1):
            
                                bs_ind[j] = thick_fil_0[f"cb_nearest_a_n[x_{j}]"][cb_ind]                                       
                                stretch[j] = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind[j]] 
                                cb_stretch_0 = kd.get_stretch_interval(stretch[j])
                                
                                complete_transition[cb_iso_0-1,idx,cb_stretch_0] += 1
                                
                        elif cb_1_type == 'S':
                            
                            # Check is cb_ind is odd -> no complete transition to account for
                            # Complete transition if cb_ind is even
                            # Check if the SRX dimers transition at the same time
                            
                            break
                            
                ### Determine all potential transitions depending on state type (S, D and A) ###
                
                # SRX (S)
                
                if cb_0_type == 'S':
                    
                    # Check is cb_ind is odd -> no potential transition to account for
                    # Potential transition if cb_ind is even
                    
                    break
                
                # Attached (A)
                
                if cb_0_type == 'A': # Get CB stretch

                    thin_ind = thick_fil_0["cb_bound_to_a_f"][cb_ind]
                    bs_ind_0 = thick_fil_0["cb_bound_to_a_n"][cb_ind]
                    thin_fil = hs_0["thin"][thin_ind]   # Correction 
                    
                    stretch = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind_0] 
                    
                    cb_stretch_0 = kd.get_stretch_interval(stretch)
                    
                    for trans in m_kinetics[cb_iso_0-1][cb_state_0-1]["transition"]: 
                        
                        # All transitions from an attached state are possible
                        
                        idx_pot = trans["index"]                                               
                        potential_transition[cb_iso_0-1, idx_pot, cb_stretch_0] += 1 
                        
                # Detached (D)
                
                elif cb_0_type == 'D':
                    
                    # Get array of nearest BS and calculate stretch array
                
                    bs_ind = np.zeros(2*adj_bs + 1, dtype = int)
                    stretch = np.zeros(2*adj_bs + 1)
                    
                    thin_ind = thick_fil_0["cb_nearest_a_f"][cb_ind]
                    thin_fil = hs_0["thin"][thin_ind]
                    
                    for j in range(0, 2*adj_bs+1):
    
                        bs_ind[j] = thick_fil_0[f"cb_nearest_a_n[x_{j}]"][cb_ind]                                       
                        stretch[j] = thick_fil_0["cb_x"][cb_ind] - thin_fil["bs_x"][bs_ind[j]] 
                   
                    for trans in m_kinetics[cb_iso_0-1][cb_state_0-1]["transition"]:
                        
                        new_state = trans["to"]
                        idx_pot = trans["index"]
                        
                        if m_kinetics[cb_iso_0-1][new_state-1]["state_type"]== 'A':
                            
                            # Check if the potential binding sites are available 
                            
                            for k in range(0, 2*adj_bs +1):
                                
                                bs_availability = thick_fil_0[f"cb_nearest_a_n_states[x_{k}]"][cb_ind] 
                            
                                if bs_availability == 2: # bs is free and available for attachment                                                                                
                                    cb_stretch_0 = kd.get_stretch_interval(stretch[k])
                                    potential_transition[cb_iso_0-1, idx_pot, cb_stretch_0] += 1 

                        elif m_kinetics[cb_iso_0-1][new_state-1]["state_type"]== 'S':
                            
                            # Check is cb_ind is odd -> no potential transition to account for
                            # Potential transition if cb_ind is even
                            break 
                        
                        else: # Transition to a "D" state is always possible
                            for k in range(0, 2*adj_bs +1):   
                                cb_stretch_0 = kd.get_stretch_interval(stretch[k])
                                potential_transition[cb_iso_0-1, idx_pot, cb_stretch_0] += 1                                                           
                
                else:
                    raise RuntimeError(f"State #{cb_state_0} is neither type A, D or S")
                                
              
        hs_0 = hs_1 

    
    for iso in range(0, len(m_kinetics)):
        for trans in range(0, max_no_of_trans):
            for stretch_bin in range(0,kd.NB_INTER-1):
                if potential_transition[iso, trans, stretch_bin] == 0:          
                    prob[iso, trans, stretch_bin] = 0
                else:
                    prob[iso, trans, stretch_bin] = complete_transition[iso, trans, stretch_bin]/potential_transition[iso, trans, stretch_bin]
                    if prob[iso, trans, stretch_bin] >= 1:
                        print("probability greater than 1")
                    else:
                        calculated_rate[iso, trans, stretch_bin] = -np.log(1.0 - prob[iso, trans, stretch_bin]) / time_step
                        
    ### Calculate 95% CI ###
    
    # p = (complete_transition + 2)/(potential_transition + 4)
    
    # W = 2 * np.sqrt( (p * (1 - p))/ (potential_transition + 4) )

    # conf_interval_pos = -np.log(1.0 - (p + W)) / time_step    
    # conf_interval_neg = -np.log(1.0 - ( p - W)) / time_step
   
    x = np.arange(kd.X_MIN,kd.X_MAX, kd.X_STEP)

    kd.calculate_rate_from_m_kinetics(m_kinetics, pd.MODEL_FILE)  
    
    for iso in range(0, len(m_kinetics)):
     
        rates_file = os.path.join(pd.OUTPUT_DIR, f"rate_equations_iso_{iso}.txt")  
        rate_values = np.loadtxt(rates_file)
        stretches = rate_values[:, 0]
        #node_force = rate_values[:,1]
        
        title = ["Attachment rate", "Detachment rate"]
        
        for i in range(0, max_no_of_trans):
    
            plt.figure()
            plt.plot(x, calculated_rate[iso, i, :], label = "calculated rate")
            plt.ylabel("Rate")
            plt.xlabel("CB stretch [nm]")
            plt.plot(stretches, rate_values[:, i + 2], label = "rate law")  
            plt.legend()
            plt.title(title[i])
            plt.show()  
            
        # for i in range(0, max_no_of_trans):
            
        #     # Get the asymmetrical y error bars from calculated confidence intervals.
        #     y_err = [
        #       abs(calculated_rate[iso, i, :] - conf_interval_neg[iso, i, :]),
        #       abs(calculated_rate[iso, i, :] + conf_interval_pos[iso, i, :])
        #     ]

        #     plt.figure()
        #     plt.plot(stretches, rate_values[:, i + 2], label = "rate law", color = "tab:orange")  
        #     plt.bar(stretches, calculated_rate[iso, i, :], yerr=y_err, label = "calculated rate")
        #     plt.ylabel("Rate")
        #     plt.xlabel("CB stretch [nm]")

        #     plt.legend()
        #     plt.show() 
    
    return calculated_rate


calc_rate = compute_m_kinetics_rate()

    
# from scipy.optimize import curve_fit

# def fit_exponential_recovery(x, y, n=1):
#     """ Fits exponential recovery with a single exponential of form y = offset + amp*(1 - exp(-k*x)) to y data """
    
#     if n==1:
#         def y_single_exp(x_data, offset, amp, k):
#             y = np.zeros(len(x_data))
#             for i,x in enumerate(x_data):
#                 y[i] = offset + amp*(1 - np.exp(-k*x))
#             return y
        
#         popt, pcov = curve_fit(y_single_exp, x, y,
#                                [y[0], y[-1]-y[0], (1/(0.2*np.amax(x)))])
        
#         d = dict()
#         d['offset'] = popt[0]
#         d['amp'] = popt[1]
#         d['k'] = popt[2]
#         d['x_fit'] = np.linspace(x[0], x[-1], 1000)
#         d['y_fit'] = y_single_exp(d['x_fit'], *popt)
        
#         return d


# results_file = os.path.join(pd.OUTPUT_DIR, "results.txt") 
# results_file = "C:/Users/sako231/OneDrive - University of Kentucky/Documents/Python Scripts/Kinetic_test_new/run_with_log/sim_output/results.txt" 
# results = np.loadtxt(results_file, skiprows=2, usecols = (0,13))

# time = results[:,0]
# m_pop_1 = results[:,1]

# plt.figure()

# plt.plot(time, m_pop_1)

# fit_data = fit_exponential_recovery(time, m_pop_1, n=1)

# plt.plot(fit_data['x_fit'], fit_data['y_fit'], color = "tab:red")
# plt.title("m_pop_1")
# plt.xlabel("Time (s)")
# plt.ylabel("Population fraction")
# print(fit_data['k'])

# plt.show()





