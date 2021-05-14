# -*- coding: utf-8 -*-
"""
Created on Thu May 13 17:59:11 2021

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

def compute_a_kinetics_rate():
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
    pCa = protocol[0][1]
    calcium = pow(10,-pCa)
    
    ### Get thin kinetic and structure parameters
   
    with open(pd.MODEL_FILE, 'r') as f:
        mod = json.load(f)
           
    k_on = mod["thin_parameters"]["a_k_on"]
    k_off = mod["thin_parameters"]["a_k_off"]
    k_coop = mod["thin_parameters"]["a_k_coop"]
    #nb_of_bs_state = mod["thin_parameters"]["a_no_of_bs_states"]
    
    a_strands_per_filament = mod["thin_structure"]["a_strands_per_filament"]
    a_regulatory_units_per_strand = mod["thin_structure"]["a_regulatory_units_per_strand"]
    a_bs_per_unit = mod["thin_structure"]["a_bs_per_unit"]

    
    ### Initialize transition arrays
    
    trans_off = np.zeros(3,dtype=int)
    trans_on = np.zeros(3,dtype=int)
    pot_off = np.zeros(3,dtype=int)
    pot_on = np.zeros(3,dtype=int)
    
    prob_on = np.zeros(3,dtype=float)
    rate_on = np.zeros(3,dtype=float)
    prob_off = np.zeros(3,dtype=float)
    rate_off = np.zeros(3,dtype=float)

               
    # HS at t
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])["hs_data"]
  
    for ind in range(1,len(hs_file)): 
        # HS at t + dt
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])["hs_data"]
      
        for i, thin_fil_0 in enumerate(hs_0["thin"]): # Loop over all thin filaments
            
            thin_fil_1 = hs_1["thin"][i] 
   
            for unit_ind, unit_state in enumerate(thin_fil_0["unit_status"]): # loop over all units
                
                unit_state_0 = thin_fil_0["unit_status"][unit_ind] # Unit at time t
                unit_state_1 = thin_fil_1["unit_status"][unit_ind] # Unit at time t + dt
                
                bs_idx_start = (unit_ind // a_strands_per_filament) * a_bs_per_unit * a_strands_per_filament + (unit_ind %2) 
                
                # Evaluate the neighbouring units status
                
                coop_left = 0
                coop_right = 0
                coop_idx = 0
                
                if unit_ind > 1:
                    coop_left += thin_fil_0["unit_status"][unit_ind - 2] - 1                                       
                    
                if unit_ind < a_strands_per_filament * a_regulatory_units_per_strand - 2:
                    coop_right += thin_fil_0["unit_status"][unit_ind + 2] - 1
                                                        
                coop_idx = coop_left + coop_right
                
                # Evaluate unit status
                                      
                if unit_state_0 == 1: # Unit is OFF                  
                    pot_on[coop_idx] += 1 # Unit can always turn ON
                    
                    if unit_state_1 == 2: # Unit activates
                        trans_on[coop_idx] +=1 
                        
                elif unit_state_0 == 2: # Unit is ON

                    # Unit may deactivate only if there is no CB bound within the unit
                    
                    unit_occupied = False
                    
                    for j in range(0,a_bs_per_unit): # Loop over all bs from the unit
                        
                        bs_idx = bs_idx_start + 2*j   
                        
                        if thin_fil_0["bound_to_m_f"][bs_idx] != -1:
                            unit_occupied = True
                            break

                    if unit_occupied == False: # Unit is freeD FROM DESIRE, MIND AND SENSES PURIFIED
                        pot_off[coop_idx] += 1
                            
                    if unit_state_1 == 1: # Unit deactivates
                        trans_off[coop_idx] += 1
                        
                else:
                    print("Unit state is neither 1 or 2")
     
        hs_0 = hs_1
        
    ### Calculate the rate 
    
    for ind in range(0,3):

        if trans_off[ind] != 0:        
            prob_off[ind] = trans_off[ind]/pot_off[ind]
        else:
            print("potential number of on-to-off transition is null")
            prob_off[ind] = 0
        
        if trans_on[ind] != 0:        
            prob_on[ind] = trans_on[ind]/pot_on[ind]
        else:
            print("potential number of off-to-on transition is null")
            prob_on[ind] = 0

        rate_on[ind] = -np.log(1.0 - prob_on[ind]) / time_step / calcium    
        rate_off[ind] = -np.log(1.0 - prob_off[ind]) / time_step 
        
    ### Compute the 95% CI 
    
    rate_conf_inter_off_neg = np.zeros(3, dtype = float)
    rate_conf_inter_off_pos = np.zeros(3, dtype = float)
    
    rate_conf_inter_on_neg = np.zeros(3, dtype = float)
    rate_conf_inter_on_pos = np.zeros(3, dtype = float)
    
    p_on = (trans_on + 2)/(pot_on + 4)
    W_on = 2*np.sqrt(p_on * (1 - p_on)/(pot_on + 4))

    conf_inter_on_neg = p_on - W_on
    conf_inter_on_pos = p_on + W_on

    rate_conf_inter_on_neg = -np.log(1.0 -conf_inter_on_neg) / time_step
    rate_conf_inter_on_pos = -np.log(1.0 -conf_inter_on_pos) / time_step

    
    p_off = (trans_off + 2)/(pot_off + 4)
    W_off = 2*np.sqrt(p_off * (1 - p_off)/(pot_off + 4))

    conf_inter_off_neg = p_off - W_off
    conf_inter_off_pos = p_off + W_off

    rate_conf_inter_off_neg = -np.log(1.0 -conf_inter_off_neg) / time_step
    rate_conf_inter_off_pos = -np.log(1.0 -conf_inter_off_pos) / time_step
        
    ### Plot the results
    
    coop = [1,2,3]
    
    k_on = [k_on, k_on * (1 + k_coop), k_on * ( 1 + 2 * k_coop)]
    k_off = [k_off * ( 1 + 2 * k_coop), k_off * (1 + k_coop), k_off]
    
    # Transition from OFF to ON
    plt.figure()
    plt.errorbar(coop, rate_on, yerr = [rate_conf_inter_on_neg, rate_conf_inter_on_pos], label = "calculated rate", 
              ecolor = "tab:grey", fmt='s', markersize = 8, markerfacecolor='tab:grey', markeredgecolor = "tab:grey", zorder=1)
    plt.plot(coop, k_on, 'ro', label = "model rate", markersize = 6, zorder=2)
    plt.ylabel("Rate (s$^{-1}$)")
    #plt.xlabel("Cooperativity level")
    plt.xlim([0,4])
    plt.xticks([1,2,3], labels = ["no coop", "coop", "double coop"])
    #plt.yticks([])
    #plt.legend(loc = "upper left")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.show()
    
    # Transition from ON to OFF
        
    plt.figure()
    plt.errorbar(coop, rate_off, yerr = [rate_conf_inter_off_neg, rate_conf_inter_off_pos], label = "calculated rate", 
                  ecolor = "tab:grey", fmt='s', markersize = 8, markerfacecolor='tab:grey', markeredgecolor = "tab:grey", zorder=1)
    plt.plot(coop, k_off, 'ro', label = "model rate", markersize = 6, zorder=2)
    plt.ylabel("Rate (s$^{-1}$)")
    #plt.xlabel("Cooperativity level")
    plt.xlim([0,4])
    plt.xticks([1,2,3], labels = ["double coop", "coop", "no coop"])
    #plt.yticks([])
    #plt.legend(loc = "upper left")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.show()
              
compute_a_kinetics_rate()
