# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:08:32 2021

@author: sako231
"""

import os, sys
# import json

ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

import path_definitions as pd

from natsort import natsorted

import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import numpy as np

def get_hs_thin_and_thick_errors(dump_folder):
    
    hs_file = []

    for filenames in os.listdir(dump_folder):
        
        filenames = os.path.join(dump_folder, filenames)
        hs_file.append(filenames)

    hs_file = natsorted(hs_file) # sorting the dump files
    
    thick_error = []
    thin_error = []

    for ind in range(1,len(hs_file)): # loop over the hs
    
        if ind%100 == 0: # too many files

            hs = half_sarcomere.half_sarcomere(hs_file[ind])["hs_data"]
          
            thick_error.append(get_thick_node_max_error(hs))
            thin_error.append(get_thin_node_max_error(hs))
            
    plt.figure()
    plt.plot(thick_error, color = "tab:red", label = "thick error")
    plt.plot(thin_error, color = "tab:blue", label = "thin error")
    plt.legend()
     

def get_thick_node_max_error(hs):
    """
    Get the maximum node error for thick filaments
    
    """
    
    hs_length = hs["hs_data"]["hs_length"]
    thick_node_err = []

    for i, thick_fil in enumerate(hs["thick"]): # Loop over all thick filaments 
    
        
        node_err = 0.0
        
        k_m = thick_fil["m_k_stiff"]
        m_rl = thick_fil["m_inter_crown_rest_length"]
        cbs_per_node = thick_fil["m_cbs_per_node"]
        
        # First node
        
        x_0 = thick_fil["cb_x"][0]
        x_1 = thick_fil["cb_x"][0 + cbs_per_node]
        
        lamb = thick_fil["m_lambda"]
        
        residual = k_m * (x_0 - x_1 - m_rl) - k_m * (hs_length - lamb - x_0 - m_rl)
        residual += get_m_cb_force(hs, thick_fil, 0, 0 + cbs_per_node)
        # residual += get_titin_force(hs, thick_fil, 0, 0 + cbs_per_node)
        # residual += get_pc_force(hs, thick_fil, 0, 0 + cbs_per_node)
        
        node_err = max(node_err, abs(residual))

        # Intermediate nodes
        
        for j in range(1, hs["hs_data"]["m_nodes_per_thick_filament"] - 1):
            
            x_0 = thick_fil["cb_x"][(j - 1) * cbs_per_node] # left node
            x_1 = thick_fil["cb_x"][j * cbs_per_node]       # center node
            x_2 = thick_fil["cb_x"][(j + 1) * cbs_per_node] # right node
            
            residual = k_m * (x_0 - x_1 - m_rl) - k_m * (x_1- x_2 - m_rl)
            residual += get_m_cb_force(hs, thick_fil, j, j + cbs_per_node)
            
            node_err = max(node_err, abs(residual))
            
        # Last node
        
        x_0 = thick_fil["cb_x"][(hs["hs_data"]["m_nodes_per_thick_filament"] - 2) * cbs_per_node] # second to last node
        x_1 = thick_fil["cb_x"][-1] # last element
        
        residual = k_m * (x_0 - x_1 - m_rl)
        residual += get_m_cb_force(hs, thick_fil, hs["hs_data"]["m_nodes_per_thick_filament"] - 2, hs["hs_data"]["m_nodes_per_thick_filament"] - 1)
            
        node_err = max(node_err, abs(residual)) 

        thick_node_err.append(node_err)
        
    max_thick_node_error = max(thick_node_err)/k_m # convert to nm
        
    return max_thick_node_error

def get_thin_node_max_error(hs):
    """
    Get the maximum node error for thin filaments
    
    """

    thin_node_err = []

    for i, thin_fil in enumerate(hs["thin"]): # Loop over all thin filaments 
    
        nearest_thicks = get_nearest_thick_fil(hs, thin_fil)
        
        node_err = 0.0
        
        k_a = thin_fil["a_k_stiff"]
        a_rl = thin_fil["a_inter_bs_rest_length"]
        bs_per_node = thin_fil["a_bs_per_node"]
        
        # First node
        
        x_0 = thin_fil["bs_x"][0]
        x_1 = thin_fil["bs_x"][0 + bs_per_node]
                
        residual = k_a * (x_1 - x_0 - a_rl)
        residual -= get_a_cb_force(hs, thin_fil, 0, 0 + bs_per_node)
        residual -= get_a_titin_force(hs, thin_fil, 0, nearest_thicks)
        # residual -= get_pc_force(hs, thick_fil, 0, 0 + cbs_per_node)
        
        node_err = max(node_err, abs(residual)) 

        # Intermediate nodes
        
        for j in range(1, hs["hs_data"]["a_nodes_per_thin_filament"] - 1):
            
            x_0 = thin_fil["bs_x"][(j - 1) * bs_per_node] # left node
            x_1 = thin_fil["bs_x"][j * bs_per_node]       # center node
            x_2 = thin_fil["bs_x"][(j + 1) * bs_per_node] # right node
            
            residual = k_a * (x_1 - x_0 - a_rl) - k_a * (x_2- x_1 - a_rl)
            residual -= get_a_cb_force(hs, thin_fil, j, j + bs_per_node)
            residual -= get_a_titin_force(hs, thin_fil, j, nearest_thicks)
            
            node_err = max(node_err, abs(residual)) 
            
        # Last node
        
        x_0 = thin_fil["bs_x"][(hs["hs_data"]["a_nodes_per_thin_filament"]  - 2) * bs_per_node] # second to last node
        x_1 = thin_fil["bs_x"][(hs["hs_data"]["a_nodes_per_thin_filament"]  - 1) * bs_per_node] # last element
        
        residual = k_a * (x_1 - x_0 - a_rl)
        residual -= get_a_cb_force(hs, thin_fil, hs["hs_data"]["a_nodes_per_thin_filament"] - 1, hs["hs_data"]["a_nodes_per_thin_filament"])
        residual -= get_a_titin_force(hs, thin_fil, hs["hs_data"]["a_nodes_per_thin_filament"]- 1, nearest_thicks)
            
        node_err = max(node_err, abs(residual))   

        thin_node_err.append(node_err)
        
    max_thin_node_error = max(thin_node_err)/k_a # convert to nm
        
    return max_thin_node_error
           

def get_m_cb_force(hs, thick_fil, idx_1, idx_2):
    """
    Calculate bound cbs force contribution to thick node forces
    """
    
    cb_force = 0.0
    
    for cb_idx in range(idx_1, idx_2):
        
        if thick_fil["cb_bound_to_a_f"][cb_idx] >= 0:
            
          cb_state = thick_fil["cb_state"][cb_idx]
          cb_extension = hs["hs_data"]["cb_extensions"][cb_state - 1]
          
          thin_fil_id = thick_fil["cb_bound_to_a_f"][cb_idx]
          thin_fil = hs["thin"][thin_fil_id]
          
          bs_idx = thick_fil["cb_bound_to_a_n"][cb_idx]
          bs_x = thin_fil["bs_x"][bs_idx]
          
          cb_x = thick_fil["cb_x"][cb_idx]
    
          delta_x = cb_x + cb_extension - bs_x
    
          cb_force += thick_fil["m_k_cb"] * delta_x
  
    return cb_force

def get_a_cb_force(hs, thin_fil, idx_1, idx_2):
    """
    Calculate bound cbs force contribution to thick node forces
    """
    cb_force = 0.0
      
    for bs_idx in range(idx_1, idx_2):
        
      if thin_fil["bound_to_m_f"][bs_idx] >= 0:
          
        thick_fil_id = thin_fil["bound_to_m_f"][bs_idx]
        thick_fil = hs["thick"][thick_fil_id]
        
        cb_idx = thin_fil["bound_to_m_n"][bs_idx]
        cb_x = thick_fil["cb_x"][cb_idx]
        
        cb_state = cb_state = thick_fil["cb_state"][cb_idx]
        cb_extension = hs["hs_data"]["cb_extensions"][cb_state - 1]
        
        bs_x = thin_fil["bs_x"][bs_idx]
        
        delta_x = cb_x + cb_extension - bs_x
  
        cb_force += thick_fil["m_k_cb"] * delta_x
    
      return cb_force
  
    
def get_a_titin_force(hs, thin_fil, bs_node, nearest_thicks, bs_per_node=2, cb_per_node=6):
      """
      Returns titin contribution to thin node force
    
      """
      titin_force = 0
      
      if bs_node == hs["titin"]["t_attach_a_node"] - 1:
          
          bs_idx = bs_node * bs_per_node 
          bs_x = thin_fil["bs_x"][bs_idx]
          
          cb_idx = cb_per_node * ( hs["titin"]["t_attach_m_node"] - 1 )
    
          for thick_fil_id in nearest_thicks:
              
            thick_fil = hs["thick"][thick_fil_id]
            
            cb_x = thick_fil["cb_x"][cb_idx]

            delta_x = cb_x - bs_x 
            
            titin_force += hs["titin"]["t_k_stiff"] * (delta_x - hs["titin"]["t_slack_length"])
    
      return titin_force

def get_m_titin_force(hs, thick_fil, cb_node, bs_per_node=2, cb_per_node=6):
      """
      Returns titin contribution to thick node force
    
      """
      titin_force = 0
      
      if cb_node == hs["titin"]["t_attach_m_node"] - 1:
          
          cb_idx = cb_node * cb_per_node 
          cb_x = thick_fil["cb_x"][cb_idx]
          
          bs_idx = bs_per_node * (hs["titin"]["t_attach_a_node"] - 1)

          for thin_id in thick_fil["nearest_actin_filaments"]:
              
            thin_fil = hs["thin"][thin_id]
            
            bs_x = thin_fil["bs_x"][bs_idx]
            
            delta_x = cb_x - bs_x 
            
            titin_force += hs["titin"]["t_k_stiff"] * (delta_x - hs["titin"]["t_slack_length"])
    
      return titin_force
  
def get_nearest_thick_fil(hs, thin_fil):
    """
    Returns array of the 3 nearest thick filaments  
    
    """
    nearest_thicks = []
    
    for thick_fil in hs["thick"]:
        
      if thin_fil["thin_id"] in thick_fil["nearest_actin_filaments"]:
          
          nearest_thicks.append(thick_fil["thick_id"])

    return nearest_thicks
    
get_hs_thin_and_thick_errors(pd.HS_STATUS_FOLDER)
                                                
