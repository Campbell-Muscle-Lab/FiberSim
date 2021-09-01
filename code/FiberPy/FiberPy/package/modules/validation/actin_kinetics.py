import os
import json

import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': "Arial"})
plt.rcParams.update({'font.size': 14}) 

# Add the half-sarcomere package

# ROOT = os.path.dirname(__file__)
# MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
# sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere


def compute_rate(model_file, protocol_file, dump_folder, output_folder):
    """Approximates the governing rate functions for actin"""

    ### Get the hs_status dump files ###

    hs_file = []

    for filenames in os.listdir(dump_folder):
        filenames = os.path.join(dump_folder, filenames)

        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # sorting the dump files
    
    ### Get the time step and calcium concentration ###

    protocol = np.loadtxt(protocol_file, skiprows=1)
    time_step = protocol[0][0]
    pCa = protocol[0][1]
    calcium = pow(10,-pCa)
    
    ### Get thin kinetic and structure parameters
   
    with open(model_file, 'r') as f:
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
    hs_0 = half_sarcomere.half_sarcomere(hs_file[0])
  
    for ind in range(1,len(hs_file)): 
        # HS at t + dt
        hs_1 = half_sarcomere.half_sarcomere(hs_file[ind])
      
        for i, thin_fil_0 in enumerate(hs_0["thin"]): # Loop over all thin filaments
            
            thin_fil_1 = hs_1["thin"][i] 
   
            for unit_ind, unit_state in enumerate(thin_fil_0["unit_status"]): # loop over all units
                
                unit_state_0 = thin_fil_0["unit_status"][unit_ind] # Unit at time t
                unit_state_1 = thin_fil_1["unit_status"][unit_ind] # Unit at time t + dt
                
                bs_idx_start = (unit_ind // a_strands_per_filament) * a_bs_per_unit * a_strands_per_filament + (unit_ind %2) 
            
                
                # Evaluate unit status
                                      
                if unit_state_0 == 1: # Unit is OFF   
                
                    # Evaluate the neighbouring units status
                
                    coop_left = 0
                    coop_right = 0
                    coop_idx = 0
                    
                    if unit_ind > 1:
                        coop_left += thin_fil_0["unit_status"][unit_ind - 2] - 1                                       
                        
                    if unit_ind < a_strands_per_filament * a_regulatory_units_per_strand - 2:
                        coop_right += thin_fil_0["unit_status"][unit_ind + 2] - 1
                                                            
                    coop_idx = coop_left + coop_right
                
                    pot_on[coop_idx] += 1 # Unit can always turn ON
                    
                    if unit_state_1 == 2: # Unit activates
                        trans_on[coop_idx] +=1 
                        
                elif unit_state_0 == 2: # Unit is ON
                
                    # Evaluate the neighbouring units status
                
                    coop_left = 0
                    coop_right = 0
                    coop_idx = 0
                    
                    if unit_ind > 1:
                        if thin_fil_0["unit_status"][unit_ind - 2] == 1:
                            coop_left = 1                                      
                        
                    if unit_ind < a_strands_per_filament * a_regulatory_units_per_strand - 2:
                        if thin_fil_0["unit_status"][unit_ind + 2] == 1:
                            coop_right = 1
                                                            
                    coop_idx = coop_left + coop_right

                    # Unit may deactivate only if there is no CB bound within the unit
                    
                    unit_occupied = 0
                    
                    for j in range(0,a_bs_per_unit): # Loop over all bs from the unit
                        
                        bs_idx = bs_idx_start + 2*j   
                        
                        if thin_fil_0["bound_to_m_f"][bs_idx] != -1:
                            unit_occupied = 1
                            break

                    if unit_occupied == 0: # Unit is freeD FROM DESIRE, MIND AND SENSES PURIFIED
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

    rate_conf_inter_on_neg = -np.log(1.0 -conf_inter_on_neg) / time_step / calcium
    rate_conf_inter_on_pos = -np.log(1.0 -conf_inter_on_pos) / time_step / calcium

    
    p_off = (trans_off + 2)/(pot_off + 4)
    W_off = 2*np.sqrt(p_off * (1 - p_off)/(pot_off + 4))

    conf_inter_off_neg = p_off - W_off
    conf_inter_off_pos = p_off + W_off

    rate_conf_inter_off_neg = -np.log(1.0 -conf_inter_off_neg) / time_step
    rate_conf_inter_off_pos = -np.log(1.0 -conf_inter_off_pos) / time_step
        
    ### Plot the results
    
    coop = [1,2,3]
    
    k_on = [k_on, k_on * (1.0 + k_coop), k_on * ( 1.0 + 2 * k_coop)]
    k_off = [k_off, k_off * (1.0 + k_coop), k_off * (1.0 + 2 * k_coop)]
    
    # Transition from OFF to ON
    
    output_image_file = os.path.join(output_folder, "inactive_to_active.png")
    
    y_err = [abs(rate_conf_inter_on_neg - rate_on), 
             abs(rate_conf_inter_on_pos - rate_on)]
        
    plt.figure(constrained_layout=True)
    plt.errorbar(coop, rate_on, yerr = y_err, label = "calculated rate", 
              ecolor = "tab:grey", fmt='s', markersize = 8, markerfacecolor='tab:grey', markeredgecolor = "tab:grey", zorder=1)
    plt.plot(coop, k_on, 'ro', label = "model rate", markersize = 6, zorder=2)
    #plt.plot(coop, rate_on, 'bs', label = "calculated rate", markersize = 6, zorder=2)
    plt.ylabel("Rate (M$^{-1}$ s$^{-1}$)")
    plt.xlim([0,4])
    plt.xticks([1,2,3], labels = ["no coop", "coop", "double coop"])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.title("Inactive to active transition rates")
    plt.legend()
    plt.savefig(output_image_file, dpi = 150)
    
    # Transition from ON to OFF
    
    output_image_file = os.path.join(output_folder, "active_to_inactive.png")
    
    y_err = [abs(rate_conf_inter_off_neg - rate_off), 
             abs(rate_conf_inter_off_pos - rate_off)]
        
    plt.figure(constrained_layout=True)
    plt.errorbar(coop, rate_off, yerr = y_err, label = "calculated rate", 
                  ecolor = "tab:grey", fmt='s', markersize = 8, markerfacecolor='tab:grey', markeredgecolor = "tab:grey", zorder=1)
    plt.plot(coop, k_off, 'ro', label = "model rate", markersize = 6, zorder=2)
    #plt.plot(coop, rate_off, 'bs', label = "calculated rate", markersize = 6, zorder=2)
    plt.ylabel("Rate (s$^{-1}$)")
    plt.xlim([0,4])
    plt.xticks([1,2,3], labels = ["no coop", "coop", "double coop"])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.title("Active to inactive transition rates")
    plt.legend()
    plt.savefig(output_image_file, dpi = 150)
              
