# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:08:32 2021

@author: sako231
"""

import os
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

from modules.half_sarcomere import half_sarcomere


def compute_force_balance(dump_precision_list, dump_folder_list, output_folder) :
    
    max_thick_node_err = []
    max_thin_node_err = []
    
    total_force_err = []
    
    # cb = np.zeros((2, 25, 54), dtype = float)
    # cb_stretch = np.zeros((2, 25, 53), dtype = float)
    # cb_stretch_2 = np.zeros((2, 25, 52), dtype = float)
    
    # bs = np.zeros((2, 50, 189), dtype = float)
    # bs_stretch = np.zeros((2, 50, 188), dtype = float)
    # bs_stretch_2 = np.zeros((2, 50, 187), dtype = float)
    
    # idx_1 = 0
    
    for dump_folder in dump_folder_list: # Loop over each dump folder
    
        hs_file = []
        
        for filenames in os.listdir(dump_folder):
            
            filenames = os.path.join(dump_folder, filenames)
            hs_file.append(filenames)
        
        hs_file = natsorted(hs_file) # sorting the dump files   
        
        hs = half_sarcomere.half_sarcomere(hs_file[-1]) # get last dump_file
        
        thick_err, thin_err = get_hs_thin_and_thick_errors(hs)
        
        max_thick_node_err.append(max(thick_err))
        
        max_thin_node_err.append(max(thin_err))

        total_force_err.append(check_total_force(hs))
        
    #     ### DIFFERENCE IN THICK NODE POSITIONS
        
    #     for i, thick_fil in enumerate(hs["thick"]): # Loop over all thick filaments 
        
    #         for j in range(54):
        
    #             cb[idx_1, i, j] = thick_fil["node_forces"][j]/20 #thick_fil["cb_x"][j * 6] # thick node positions
                
    #         for j in range(53):
                                
    #             cb_stretch[idx_1, i, j] = thick_fil["cb_x"][j * 6] - thick_fil["cb_x"][(j+1) * 6] # delta x
        
    #         for j in range(52):
                
    #             cb_stretch_2[idx_1, i, j]  = cb_stretch[idx_1, i, j+1]  -  cb_stretch[idx_1, i, j] 
        
    #     for i, thin_fil in enumerate(hs["thin"]): # Loop over all thin filaments 
        
    #         for j in range(189):
        
    #             bs[idx_1, i, j] = thin_fil["bs_x"][j * 2] # thin node positions
        
    #         for j in range(188):
                                
    #             bs_stretch[idx_1, i, j] = thin_fil["bs_x"][(j+1) * 2] - thin_fil["bs_x"][j * 2] # delta x
    #         for j in range(187):
                                
    #             bs_stretch_2[idx_1, i, j] = bs_stretch[idx_1, i, j+1] - bs_stretch[idx_1, i, j]
                
    #     idx_1 += 1

    # plt.figure()
    
    # for j in range(9):
        
    #     plt.plot(cb[0, j, :] - cb[1, j, :])
        
    #     plt.title("Difference in thick nodes positions between low and high tolerance dumped values")
            
    # plt.show()
    
    # for j in range(9):
        
    #     plt.plot(cb_stretch[0, j, :], 'k')
    #     plt.plot(cb_stretch[1, j, :], 'b')
        
    #     plt.title("Difference in thick springs stretches between low and high tolerance dumped values")
            
    # plt.show()
    
    # for j in range(18):
        
    #     plt.plot(bs[0, j, :] - bs[1, j, :])
        
    #     plt.title("Difference in thin nodes positions between low and high tolerance dumped values")
            
    # plt.show()
    
    # for j in range(18):
        
    #     plt.plot(bs_stretch[0, j, :], 'k')
    #     plt.plot(bs_stretch[1, j, :], 'b')
        
    #     plt.title("Difference in thin springs stretches between low and high tolerance dumped values")
            
    # plt.show()
    
    # PLOTS    
    
    output_image_file = os.path.join(output_folder, "Force_error.png")
    
    ## Check that pc nearest actin distribution is uniform
    
    # counter = np.zeros(len(hs["thin"]))
    
    # for i, thick_fil in enumerate(hs["thick"]): # Loop over all thick filaments 
    
    #     nearest = thick_fil["nearest_actin_filaments"]
        
    #     for thin_fil in nearest:
            
    #         for nearest_thin in thick_fil["pc_nearest_a_f"]:
            
    #             if thin_fil in thick_fil["pc_nearest_a_f"]:
                    
    #                 counter[thin_fil] += 1
        
    # print("NEAREST ACTIN FILAMENT FOR PC =")
    # print(counter)
    
    plt.figure()
    plt.plot(dump_precision_list, total_force_err, color = "tab:red")
    plt.xlabel("dump precision (# digits)")
    plt.xticks(dump_precision_list)
    plt.ylabel("Total force error (%)")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    #plt.yscale("log")
    
    #plt.savefig(output_image_file, dpi = 150)

    output_image_file = os.path.join(output_folder, "Filaments_error.png")
    
    plt.figure()
    plt.plot(dump_precision_list, max_thick_node_err, color = "tab:red", label = "average max thick node error")
    plt.plot(dump_precision_list, max_thin_node_err, color = "tab:blue", label = "average max thin node error")
    plt.xlabel("dump precision (# digits)")
    plt.xticks(dump_precision_list)
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel("Error [nm]")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()
    
    plt.savefig(output_image_file, dpi = 150)

def get_hs_thin_and_thick_errors(hs, show_error = False):
     
    thick_error = get_thick_node_error(hs)
    thin_error = get_thin_node_error(hs)
            
    if show_error:
        plt.figure()
        plt.plot(thick_error, color = "tab:red", label = "thick error")
        plt.plot(thin_error, color = "tab:blue", label = "thin error")
        plt.xlabel("# fil")
        plt.ylabel(" error (nm)")
        plt.yscale("log")
        plt.legend()
        plt.show()
    return thick_error, thin_error
     

def get_thick_node_error(hs):
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
        residual += get_m_pc_force(hs, thick_fil, 0)
        residual += get_m_titin_force(hs, thick_fil, 0)
        
        node_err = max(node_err, abs(residual))

        # Intermediate nodes
        
        for j in range(1, hs["hs_data"]["m_nodes_per_thick_filament"] - 1):
            
            x_0 = thick_fil["cb_x"][(j - 1) * cbs_per_node] # left node
            x_1 = thick_fil["cb_x"][j * cbs_per_node]       # center node
            x_2 = thick_fil["cb_x"][(j + 1) * cbs_per_node] # right node
            
            residual = -k_m * (x_0 - x_1 - m_rl) + k_m * (x_1- x_2 - m_rl) # TO CHECK
            residual += get_m_cb_force(hs, thick_fil, j * cbs_per_node, (j + 1) * cbs_per_node)
            residual += get_m_pc_force(hs, thick_fil, j)
            residual += get_m_titin_force(hs, thick_fil, j) 
            
            # node_force = k_m * (x_0 - x_1 - m_rl) 
            # real_node = thick_fil["node_forces"][j]
            # diff = node_force - real_node
            # print("node_force_diff =", diff)

            
            node_err = max(node_err, abs(residual))
            
        # Last node
        
        x_0 = thick_fil["cb_x"][(hs["hs_data"]["m_nodes_per_thick_filament"] - 2) * cbs_per_node] # second to last node
        x_1 = thick_fil["cb_x"][-1] # last element
        
        residual = - k_m * (x_0 - x_1 - m_rl) # FIXED ERROR
        residual += get_m_cb_force(hs, thick_fil, (hs["hs_data"]["m_nodes_per_thick_filament"] - 1) * cbs_per_node, (hs["hs_data"]["m_nodes_per_thick_filament"] - 0) * cbs_per_node)
        residual += get_m_pc_force(hs, thick_fil, hs["hs_data"]["m_nodes_per_thick_filament"] - 1)
        residual += get_m_titin_force(hs, thick_fil, hs["hs_data"]["m_nodes_per_thick_filament"] - 1) 
        
        node_err = max(node_err, abs(residual)) 
        
        thick_node_err.append(node_err)
        
    thick_node_err = [x/k_m for x in thick_node_err] # convert to nm
        
    return thick_node_err

def get_thin_node_error(hs):
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
                
        residual = k_a * (x_0 - a_rl) - k_a * (x_1 - x_0 - a_rl) # FIXED ERROR
        residual -= get_a_cb_force(hs, thin_fil, 0, 0 + bs_per_node)
        residual -= get_a_pc_force(hs, thin_fil, 0, 0 + bs_per_node)
        residual -= get_a_titin_force(hs, thin_fil, 0, nearest_thicks)
        
        node_err = max(node_err, abs(residual)) 
        

        # Intermediate nodes
        
        for j in range(1, hs["hs_data"]["a_nodes_per_thin_filament"] - 1):
            
            x_0 = thin_fil["bs_x"][(j - 1) * bs_per_node] # left node
            x_1 = thin_fil["bs_x"][j * bs_per_node]       # center node
            x_2 = thin_fil["bs_x"][(j + 1) * bs_per_node] # right node
                        
            
            residual = k_a * (x_1 - x_0 - a_rl) - k_a * (x_2- x_1 - a_rl)
            residual -= get_a_cb_force(hs, thin_fil, j * bs_per_node, (j + 1 ) * bs_per_node)
            residual -= get_a_pc_force(hs, thin_fil, j * bs_per_node, (j + 1 ) * bs_per_node)
            residual -= get_a_titin_force(hs, thin_fil, j, nearest_thicks)
            
            node_err = max(node_err, abs(residual)) 
            
        # Last node
        
        x_0 = thin_fil["bs_x"][(hs["hs_data"]["a_nodes_per_thin_filament"]  - 2) * bs_per_node] # second to last node
        x_1 = thin_fil["bs_x"][(hs["hs_data"]["a_nodes_per_thin_filament"]  - 1) * bs_per_node] # last element
        
        residual = k_a * (x_1 - x_0 - a_rl)
        residual -= get_a_cb_force(hs, thin_fil, (hs["hs_data"]["a_nodes_per_thin_filament"] - 1) * bs_per_node, hs["hs_data"]["a_nodes_per_thin_filament"] * bs_per_node)
        residual -= get_a_pc_force(hs, thin_fil, (hs["hs_data"]["a_nodes_per_thin_filament"] - 1) * bs_per_node, hs["hs_data"]["a_nodes_per_thin_filament"] * bs_per_node)
        residual -= get_a_titin_force(hs, thin_fil, hs["hs_data"]["a_nodes_per_thin_filament"]- 1, nearest_thicks)
            
        node_err = max(node_err, abs(residual))
                
        thin_node_err.append(node_err)  # store the max node error for each filament
        
    thin_node_err = [x/k_a for x in thin_node_err] # convert to nm
        
    return thin_node_err
           

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
          
        if thin_fil["bound_to_m_type"][bs_idx] == 1: # a cb is bound FIXED ERROR
        
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
  
def get_m_pc_force(hs, thick_fil, thick_node):
    """
    Calculate bound pcs force contribution to thick node forces
    """
    
    pc_force = 0.0
    
    for i, pc in enumerate(thick_fil["pc_node_index"]):
        
        if pc == thick_node: # check if there is a pc on that node
        
            if thick_fil['pc_bound_to_a_f'][i] >= 0: # check if pc is attached to actin:
                
                pc_idx = thick_node * thick_fil["m_cbs_per_node"] # get pc idx
                pc_x = thick_fil["cb_x"][pc_idx]
                
                thin_fil_id = thick_fil['pc_bound_to_a_f'][i]
                thin_fil = hs["thin"][thin_fil_id]
                
                bs_idx = thick_fil['pc_bound_to_a_n'][i]
                bs_x = thin_fil["bs_x"][bs_idx]
                
                delta_x = pc_x - bs_x
                
                pc_force += thick_fil["c_k_stiff"] * delta_x
  
    return pc_force

def get_a_pc_force(hs, thin_fil, idx_1, idx_2):
    """
    Calculate bound pcs force contribution to thin node forces
    """
    pc_force = 0.0
      
    for bs_idx in range(idx_1, idx_2):
        
        if thin_fil["bound_to_m_type"][bs_idx] == 2: # a pc is bound        
          
            thick_fil_id = thin_fil["bound_to_m_f"][bs_idx]
            thick_fil = hs["thick"][thick_fil_id]
            
            pc_idx = thick_fil["pc_bound_to_a_n"].index(bs_idx) 
            pc_node_index = thick_fil["pc_node_index"][pc_idx]
            pc_x = thick_fil["cb_x"][pc_node_index * thick_fil["m_cbs_per_node"]]
            
            #pc_idx = thin_fil["bound_to_m_n"][bs_idx]            
            #pc_x = thick_fil["cb_x"][pc_idx]
        
            bs_x = thin_fil["bs_x"][bs_idx]
        
            delta_x = pc_x - bs_x

            pc_force += thick_fil["c_k_stiff"] * delta_x

    return pc_force


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
  
def get_nearest_thick_fil(hs, thin_fil):
    """
    Returns array of the 3 nearest thick filaments  
    
    """
    nearest_thicks = []
    
    for thick_fil in hs["thick"]:
        
      if thin_fil["thin_id"] in thick_fil["nearest_actin_filaments"]:
          
          nearest_thicks.append(thick_fil["thick_id"])

    return nearest_thicks

def check_total_force(hs):
    """
    Compare the half-sarcomere force to the 
    re-calculated total force from hs dump file
    """    
    hs_force = hs["hs_data"]["hs_force"]
    hs_length = hs["hs_data"]["hs_length"]
    
    prop_fibrosis = hs["hs_data"]["prop_fibrosis"]
    prop_myofilaments = hs["hs_data"]["prop_myofilaments"]
    m_filament_density = hs["hs_data"]["m_filament_density"]
       
    computed_force = 0.0
    
    for thick_fil in hs["thick"]:
        
        k_m = thick_fil["m_k_stiff"]       
        m_rl = thick_fil["m_inter_crown_rest_length"]          
        lamb = thick_fil["m_lambda"]
        
        cb_x_m_line = thick_fil["cb_x"][0]
        
        computed_force += k_m * (hs_length - lamb - m_rl - cb_x_m_line)
        
    computed_force /= len(hs["thick"])     # average over all thick filaments
    computed_force *= (1 - prop_fibrosis) * prop_myofilaments * m_filament_density * 1e-9 # force in N/m²
       
    err_force = (hs_force - computed_force)/hs_force * 100
    
    return abs(err_force)
                                                
