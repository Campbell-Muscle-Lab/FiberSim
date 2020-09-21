# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:01:49 2020

@author: Sarah
"""
import os,sys

ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

#import sys
#sys.path.append("..") 

import numpy as np

from modules.half_sarcomere import half_sarcomere

   
def Get_error_from_hs_dump_files(file_path):
    """Gets a hs_dump_file from the dump folder, calculate the force imbalance for the thick and thin filament nodes, 
    and return the result.
  
    Parameters
    ----------
    file_path : str
        Path to the dump file.
    """    

    half = half_sarcomere.half_sarcomere(file_path)
    
    error_thick = Force_Balance_Thick(half)
    
    error_thin = Force_Balance_Thin(half)
        
    #print("Maximal error within the thick filaments = %.5f" % max(error_thick))
    
    #print("Maximal error within the thin filaments = %.5f" % max(error_thin)) 
    
    return(max(error_thin),max(error_thick))
       

def Force_Balance_Thick(half):
    """Calculate the force imbalance for the thick filaments.
  
    Parameters
    ----------
    half : str
        The half-sarcomere that is evaluated.
      
    Returns
    -------
    error_thick : list
        A list with the maximal error associated with each thick filament.
    """
    
    error_thick = [] # list with the maximal error node for each thick filament
    
    for i, thick_fil in enumerate(half.thick_fil):
        
        #print(f"Calculating the max error node for thick filament n°{i}")
        
        res = [] # list containing the nodes errors for each thick filament
        err = 0.0
                
        k_m = thick_fil.m_k_stiff 
    
        m_rl =  thick_fil.m_inter_crown_rest_length
    
        m_lambda = thick_fil.m_lambda
  
        hs_length = half.hs_length
        
        m_pos = thick_fil.cb_x
    
        no_of_cb_per_node =  thick_fil.m_cbs_per_node
        
        m_no_of_cbs = no_of_cb_per_node*half.m_nodes_per_thick_filament
        
        # First node
        
        index = 0
        node_1_index = index
        node_2_index = (index + 1) * no_of_cb_per_node
        
        err += -2*k_m*m_pos[node_1_index] + k_m*m_pos[node_2_index] + k_m*(hs_length-m_lambda)
        err-= Force_due_to_titin_binding_to_thick(thick_fil, half, index)
        err-= Force_due_to_cb_binding_to_thick(thick_fil, half, node_1_index, node_2_index)
        
        if index in thick_fil.pc_node_index:
            err-= Force_due_to_pc_binding_to_thick(thick_fil, half, index)
            
        err /= k_m
        
        res.append(err)
        
        # Intermediate nodes
    
        for index in range(1, half.m_nodes_per_thick_filament-1):
            
            err = 0.0
    
            node_0_index = (index-1) * no_of_cb_per_node
            node_1_index = index * no_of_cb_per_node 
            node_2_index = (index+1) * no_of_cb_per_node
        
            err += k_m*m_pos[node_0_index] - 2*k_m*m_pos[node_1_index] + k_m*m_pos[node_2_index]
            err -= Force_due_to_titin_binding_to_thick(thick_fil, half, index)
            err -= Force_due_to_cb_binding_to_thick(thick_fil, half, node_1_index, node_2_index)
            if index in thick_fil.pc_node_index:
                err-= Force_due_to_pc_binding_to_thick(thick_fil, half, index)
            
            err /= k_m
            
            res.append(err)
            
        # Last node
        index += 1
        
        err = 0.0
        
        node_0_index = (index-1) * no_of_cb_per_node
        node_1_index = index * no_of_cb_per_node 
        node_2_index = (index+1) * no_of_cb_per_node
            
        err += -k_m*(m_pos[node_0_index]-m_pos[node_1_index]-m_rl)
        err -= Force_due_to_titin_binding_to_thick(thick_fil, half, index)
        err -= Force_due_to_cb_binding_to_thick(thick_fil, half, node_1_index, node_2_index)
        if index in thick_fil.pc_node_index:
                err-= Force_due_to_pc_binding_to_thick(thick_fil, half, index)  
        err /= k_m
        
        res.append(err)
        
        res_abs =  [abs(ele) for ele in res] 
        
        max_error_for_this_thick = max(res_abs) # gets the maximal error node for this thick filament
        
        error_thick.append(max_error_for_this_thick) 
        
    return error_thick

def Force_Balance_Thin(half):
    """Calculate the force imbalance for the thin filaments.
  
    Parameters
    ----------
    half : 
        The half-sarcomere that is evaluated.
      
    Returns
    -------
    error_thin : list
        A list with the maximal error associated with each thin filament.
    """
    
    error_thin = [] # list with the maximal error node for each thin filament
    
    for i, thin_fil in enumerate(half.thin_fil):
        
        #print(f"Calculating the max error node for thin filament n°{i}")
        
        res = [] # list containing the nodes errors for each thin filament
        
        err = 0.0
        
        k_a = thin_fil.a_k_stiff
    
        a_rl =  thin_fil.a_inter_bs_rest_length
        
        a_pos = thin_fil.bs_x
    
        no_of_bs_per_node =  int(thin_fil.a_no_of_bs/half.a_nodes_per_thin_filament)

        
        # First node
    
        index = 0
        node_0_index = index
        node_1_index = (index + 1) * no_of_bs_per_node
               
        err += -2*k_a*a_pos[node_0_index] + k_a*a_pos[node_1_index] 
        
        err+= Force_due_to_titin_binding_to_thin(thin_fil, half, index)
        err+= Force_due_to_cb_binding_to_thin(thin_fil, half, node_0_index, node_1_index)
        err+= Force_due_to_pc_binding_to_thin(thin_fil, half, node_0_index, node_1_index)
        
        err /= k_a
        
        res.append(err)

        
        # Intermediate nodes
    
        for index in range(1, half.a_nodes_per_thin_filament-1):
            
            err = 0.0
    
            node_0_index = (index-1) * no_of_bs_per_node
            node_1_index = index * no_of_bs_per_node 
            node_2_index = (index+1) * no_of_bs_per_node
        
            err += k_a*a_pos[node_0_index] - 2*k_a*a_pos[node_1_index] + k_a*a_pos[node_2_index]
            
            err+= Force_due_to_titin_binding_to_thin(thin_fil, half, index)
            err+= Force_due_to_cb_binding_to_thin(thin_fil, half, node_1_index, node_2_index)
            err+= Force_due_to_pc_binding_to_thin(thin_fil, half, node_1_index, node_2_index)
            
            err /= k_a
            
            res.append(err)
            
        # Last node
            
        index += 1
        err = 0.0
        
        node_0_index = (index-1) * no_of_bs_per_node
        node_1_index = index * no_of_bs_per_node
        node_2_index = (index+1) * no_of_bs_per_node
            
        err += -k_a*(-a_pos[node_0_index]+a_pos[node_1_index]-a_rl)
        err += Force_due_to_titin_binding_to_thin(thin_fil, half, index)
        err+= Force_due_to_cb_binding_to_thin(thin_fil, half, node_1_index, node_2_index)
        err+= Force_due_to_pc_binding_to_thin(thin_fil, half, node_0_index, node_1_index)
        
        err /= k_a
        
        res.append(err)
        
        res_abs =  [abs(ele) for ele in res] 
        
        max_error_for_this_thin= max(res_abs) # gets the maximal error node for this thin filament
        
        error_thin.append(max_error_for_this_thin)
            
    return error_thin


def Force_due_to_titin_binding_to_thin(thin_fil, half, node_to_check):
    """Calculate the force due to titin on a thin filament node.
  
    Parameters
    ----------
    thin_fil : 
        The thin filament that is evaluated.
    half : 
        The half-sarcomere that is evaluated.
    node_to_check : int
       Index of the node that is checked. 
      
    Returns
    -------
    titin_force : float
        Titin contribution to the force balance on the node.
    """
    
    titin_force = 0.0

    
    if node_to_check + 1 == half.t_attach_a_node:
                
        for thick_fil in half.thick_fil:
             
              if thin_fil.thin_id in thick_fil.nearest_actin_filaments:
                                        
                      cb_index = thick_fil.m_cbs_per_node * half.t_attach_m_node - 1
                      m_x = thick_fil.cb_x[cb_index]
                      
                      bs_index = int(thin_fil.a_no_of_bs/half.a_nodes_per_thin_filament) * half.t_attach_a_node - 1
                      a_x = thin_fil.bs_x[bs_index]
                      
                      titin_force += half.t_k_stiff * (m_x - a_x - half.t_slack_length)
                      
                      
    return titin_force

def Force_due_to_cb_binding_to_thin(thin_fil, half, node_index_to_check_1, node_index_to_check_2):
    """Calculate the force due to cb binding on a thin filament node.
      
    Parameters
    ----------
    thin_fil : 
        The thin filament that is evaluated.
    half : 
        The half-sarcomere that is evaluated.
    node_index_to_check_1 : int
        Index of the node that is checked.
    node_index_to_check_2 : int
        Index of the node where checking stops.
      
    Returns
    -------
    cb_force : float
        cb contribution to the force balance on the node.
    """

    cb_force = 0.0
    
    for ind in range(node_index_to_check_1, node_index_to_check_2):
  
        if thin_fil.bound_to_m_n[ind] > 0:
        
            thick_fil_index = thin_fil.bound_to_m_f[ind]
            thick_fil = half.thick_fil[thick_fil_index]
        
            cb_index = thin_fil.bound_to_m_n[ind]
            m_x = thick_fil.cb_x[cb_index]
            cb_extension = half.cb_extensions[thick_fil.cb_state[cb_index]-1]
               
            a_x = thin_fil.bs_x[ind]
        
        
            cb_force += thick_fil.m_k_cb * (m_x - a_x + cb_extension)
        
    return cb_force


def Force_due_to_pc_binding_to_thin(thin_fil, half, node_index_to_check_1, node_index_to_check_2):
    """Calculate the force due to protein-C binding on a thin filament node.
  
    Parameters
    ----------
    thin_fil : 
        The thin filament that is evaluated.
    half : 
        The half-sarcomere that is evaluated.
    node_index_to_check_1 : int
        Index of the node that is checked.
    node_index_to_check_2 : int
        Index of the node where checking stops.
      
    Returns
    -------
    pc_force : float
        Protein-C contribution to the force balance on the node.
    """
        
    pc_force = 0.0
    
    for ind in range(node_index_to_check_1, node_index_to_check_2):
    
        if thin_fil.bound_to_m_f[ind] > 0: # Check for cb link
            thick_fil_id = thin_fil.bound_to_m_f[ind] # Get the thick id
            thick_fil = half.thick_fil[thick_fil_id] 
            thick_fil_node = thin_fil.bound_to_m_n[ind]//thick_fil.m_cbs_per_node # Get the myosin node number
        
            if thick_fil_node in thick_fil.pc_node_index:
                if thick_fil.pc_bound_to_a_f[thick_fil_node] > 0: # Check for pc link
                
                    a_x = thin_fil.bs_x[ind] 
                    cb_ind = thin_fil.bound_to_m_n[ind]
                    m_x = thick_fil.cb_x[cb_ind]
                    
                    pc_force += thick_fil.c_k_stiff * (m_x - a_x)
           
            
    return pc_force 




def Force_due_to_titin_binding_to_thick(thick_fil, half, node_to_check):
    """Calculate the force due to titin on a thick filament node.
  
    Parameters
    ----------
    thick_fil : 
        The thick filament that is evaluated.
    half : 
        The half-sarcomere that is evaluated.
    node_to_check : int
        Node that is checked. 
      
    Returns
    -------
    titin_force : float
        Titin contribution to the force balance on the node.
    """
    titin_force = 0.0
       
    if node_to_check + 1 == half.t_attach_m_node:
        
        for thin_fil in half.thin_fil:           
             
            if thin_fil.thin_id in thick_fil.nearest_actin_filaments:
                                        
                cb_index = thick_fil.m_cbs_per_node * half.t_attach_m_node - 1
                m_x = thick_fil.cb_x[cb_index]
                
                no_of_bs_per_node =  int(thin_fil.a_no_of_bs/half.a_nodes_per_thin_filament)
                      
                bs_index = no_of_bs_per_node * half.t_attach_a_node - 1
                a_x = thin_fil.bs_x[bs_index]
                      
                titin_force += half.t_k_stiff * (m_x - a_x - half.t_slack_length)
                      
    return titin_force

def Force_due_to_cb_binding_to_thick(thick_fil, half, node_index_to_check_1, node_index_to_check_2):
    """Calculate the force due to cb binding on a thick filament node.
  
      Parameters
      ----------
      thick_fil : 
          The thick filament that is evaluated.
      half : 
          The half-sarcomere that is evaluated.
      node_index_to_check_1 : int
         Index of the first node that is being checked 
      node_index_to_check_2 : int
         Index of the node where checking stops
        
      Returns
      -------
      cb_force : float
          cb contribution to the force balance on the node.
      """
    cb_force = 0.0
    
    for ind in range(node_index_to_check_1, node_index_to_check_2):
  
        if thick_fil.cb_bound_to_a_f[ind] > 0:
               
            thin_fil_index = thick_fil.cb_bound_to_a_f[ind]
            bs_index = thick_fil.cb_bound_to_a_n[ind]
            thin_fil = half.thin_fil[thin_fil_index]        
            a_x = thin_fil.bs_x[bs_index]
                
            m_x = thick_fil.cb_x[ind]
            cb_extension = half.cb_extensions[thick_fil.cb_state[ind]-1]
               
            cb_force += thick_fil.m_k_cb * (m_x - a_x + cb_extension)
        
    return cb_force

def Force_due_to_pc_binding_to_thick(thick_fil, half, index):
    """Calculate the force due to protein-C binding on a thick filament node.
  
    Parameters
    ----------
    thick_fil : 
        The thin filament that is evaluated.
    half : 
        The half-sarcomere that is evaluated.
    index : int
        Index of the node that is checked.
      
    Returns
    -------
    pc_force : float
        Protein-C contribution to the force balance on the node.
    """        
    pc_force = 0.0
    
    for i,pc_node in enumerate(thick_fil.pc_node_index):
        if index == pc_node and thick_fil.pc_bound_to_a_f[i] > 0:
            
            thin_fil_id = thick_fil.pc_bound_to_a_f[i]
            bs_ind = thick_fil.pc_bound_to_a_n[i]
            
            thin_fil = half.thin_fil[thin_fil_id]
            a_x = thin_fil.bs_x[bs_ind]
            
            m_x = thick_fil.cb_x[index * thick_fil.m_cbs_per_node]
            
            pc_force += thick_fil.c_k_stiff * (m_x - a_x)
            
            
    return pc_force

def Get_error_force(file_path):
    
    half = half_sarcomere.half_sarcomere(file_path)
    
    m_filament_density = 0.25e15
    
    force = half.hs_force
    computed_force = 0.0
    
    for thick_fil in half.thick_fil:
    
        k_m = thick_fil.m_k_stiff 
    
        m_rl =  thick_fil.m_inter_crown_rest_length
    
        m_lambda = thick_fil.m_lambda
  
        hs_length = half.hs_length
        
        m_pos = thick_fil.cb_x
       
        computed_force += k_m * (hs_length - m_lambda - m_rl - m_pos[0])
        
    computed_force /= len(half.thick_fil)     
    computed_force *= m_filament_density * 1e-9
    
    err_force = force-computed_force
    
    return np.abs(err_force)
            
