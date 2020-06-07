# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:11:28 2020

@author: kscamp3
"""

import numpy as np

class thin_filament():
    """Thin filament class"""
    
    def __init__(self, thin_data):
        
        self.a_no_of_bs = thin_data['a_no_of_bs']
        self.a_k_stiff = thin_data['a_k_stiff']
        self.a_inter_bs_rest_length = thin_data['a_inter_bs_rest_length']

        self.bs_x = np.array(thin_data['bs_x'])
        self.bs_angle = np.array(thin_data['bs_angle'])
        self.bs_unit = np.array(thin_data['bs_unit']).astype(int)
        self.bs_state = np.array(thin_data['bs_state']).astype(int)
        self.bs_isoform = np.array(thin_data['bs_isoform']).astype(int)
        self.bound_to_m_f = np.array(thin_data['bound_to_m_f']).astype(int)
        self.bound_to_m_n = np.array(thin_data['bound_to_m_n']).astype(int)
        self.nearest_m_f = np.array(thin_data['nearest_m_f']).astype(int)
        self.nearest_m_n = np.array(thin_data['nearest_m_n']).astype(int)

