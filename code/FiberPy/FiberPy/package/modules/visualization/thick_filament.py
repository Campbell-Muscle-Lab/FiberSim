# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:11:28 2020

@author: kscamp3
"""

import numpy as np

class thick_filament():
    """Thick filament class"""
    
    def __init__(self, thick_data):
        
        self.thick_id = thick_data['thick_id']
        self.m_no_of_cbs = thick_data['m_no_of_cbs']
        self.m_k_stiff = thick_data['m_k_stiff']
        self.m_k_cb = thick_data['m_k_cb']
        self.c_k_stiff = thick_data['c_k_stiff']
        self.m_inter_crown_rest_length = thick_data['m_inter_crown_rest_length']
        self.m_cbs_per_node = thick_data['m_cbs_per_node']
        self.m_lambda = thick_data['m_lambda']
        self.c_no_of_pcs = thick_data['c_no_of_pcs']
                
        self.m_y = thick_data['m_y']
        self.m_z = thick_data['m_z']
        
        self.nearest_actin_filaments = \
            np.array(thick_data['nearest_actin_filaments'])
        
        self.cb_x = np.array(thick_data['cb_x'])
        self.cb_angle = np.array(thick_data['cb_angle'])
        self.cb_nearest_bs_angle_diff = np.array(thick_data['cb_nearest_bs_angle_diff'])   
        self.node_forces = np.array(thick_data['node_forces'])
        self.cb_state = np.array(thick_data['cb_state']).astype(int)
        self.cb_iso = np.array(thick_data['cb_iso']).astype(int)
        self.cb_bound_to_a_f = \
            np.array(thick_data['cb_bound_to_a_f']).astype(int)
        self.cb_bound_to_a_n = \
            np.array(thick_data['cb_bound_to_a_n']).astype(int)
        self.cb_nearest_a_f = \
            np.array(thick_data['cb_nearest_a_f']).astype(int)
        self.cb_nearest_a_n = \
            np.array(thick_data['cb_nearest_a_n']).astype(int)
        self.pc_node_index = \
            np.array(thick_data['pc_node_index']).astype(int)
        self.pc_angle = np.array(thick_data['pc_angle'])
        self.pc_state = \
            np.array(thick_data['pc_state']).astype(int)
        self.pc_iso = \
            np.array(thick_data['pc_iso']).astype(int)
        self.pc_bound_to_a_f = \
            np.array(thick_data['pc_bound_to_a_f']).astype(int)
        self.pc_bound_to_a_n = \
            np.array(thick_data['pc_bound_to_a_n']).astype(int)
        self.pc_nearest_a_f = \
            np.array(thick_data['pc_nearest_a_f']).astype(int)
        self.pc_nearest_a_n = \
            np.array(thick_data['pc_nearest_a_n']).astype(int)

        self.pc_x = np.zeros((self.c_no_of_pcs,1))
        for i, ni in enumerate(self.pc_node_index):
            cb_index = ni*self.m_cbs_per_node
            self.pc_x[i] = self.cb_x[cb_index]

