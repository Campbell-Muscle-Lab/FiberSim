# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:11:28 2020

@author: kscamp3
"""

import numpy as np

class thin_filament():
    """Thin filament class"""
    
    def __init__(self, thin_data):
        
        self.thin_id = thin_data['thin_id']
        self.a_no_of_bs = thin_data['a_no_of_bs']
        self.a_bs_per_node = thin_data['a_bs_per_node']
        self.a_k_stiff = thin_data['a_k_stiff']
        self.a_inter_bs_rest_length = thin_data['a_inter_bs_rest_length']
        self.a_y = thin_data['a_y']
        self.a_z = thin_data['a_z']

        self.bs_x = np.array(thin_data['bs_x'])
        self.bs_angle = np.array(thin_data['bs_angle'])
        self.bs_unit = np.array(thin_data['bs_unit']).astype(int)
        self.bs_state = np.array(thin_data['bs_state']).astype(int)
        self.bs_isoform = np.array(thin_data['bs_isoform']).astype(int)
        self.bound_to_m_f = np.array(thin_data['bound_to_m_f']).astype(int)
        self.bound_to_m_n = np.array(thin_data['bound_to_m_n']).astype(int)
        self.bound_to_m_type = np.array(thin_data['bound_to_m_type']).astype(int)
        #self.nearest_m_f = np.array(thin_data['nearest_m_f']).astype(int)
        #self.nearest_m_n = np.array(thin_data['nearest_m_n']).astype(int)
        self.unit_status = np.array(thin_data['unit_status']).astype(int)

