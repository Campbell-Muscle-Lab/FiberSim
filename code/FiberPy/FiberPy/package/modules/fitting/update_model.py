# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 17:43:26 2023

@author: Campbell
"""

import copy

import numpy as np

def update_model(p, base_model, opt_dict, condition_number):
    
    new_model = copy.deepcopy(base_model)
    
    for (i,p_struct) in enumerate(opt_dict['parameters']):
        
        if ('condition' in p_struct):
            if not (condition_number in np.asarray(p_struct['condition'])):
                # Break out of loop if condition is specified and this
                # is not in the list
                continue
        
        param_value = return_param_value(p[i], p_struct)
        
        if (p_struct['class'] == 'lattice_parameters'):
            if (p_struct['key'] == 'viscosity'):
                new_model['lattice_parameters']['viscosity'] = param_value
            
        if (p_struct['class'].startswith('m_kinetics')):
            isotype_ind = p_struct['isotype'] - 1
            # Split key by '_'
            bits = p_struct['key'].split('_')
            state_ind = int(bits[1])-1
            trans_ind = int(bits[3])-1
            param_ind = int(bits[-1])-1
            
            y = np.asarray(new_model['m_kinetics'][isotype_ind]['scheme'][state_ind]['transition'][trans_ind]['rate_parameters'],
                           dtype=float)
            y[param_ind] = param_value
            new_model['m_kinetics'][isotype_ind]['scheme'][state_ind]['transition'][trans_ind]['rate_parameters'] = \
                y.tolist()
                
        if (p_struct['class'].startswith('c_kinetics')):
            isotype_ind = p_struct['isotype'] - 1
            # Split key by '_'
            bits = p_struct['key'].split('_')
            state_ind = int(bits[1])-1
            trans_ind = int(bits[3])-1
            param_ind = int(bits[-1])-1
            
            y = np.asarray(new_model['c_kinetics'][isotype_ind]['scheme'][state_ind]['transition'][trans_ind]['rate_parameters'],
                           dtype=float)
            y[param_ind] = param_value
            new_model['c_kinetics'][isotype_ind]['scheme'][state_ind]['transition'][trans_ind]['rate_parameters'] = \
                y.tolist()                
                
    return new_model

def return_param_value(p, p_struct):
    
    x = (p%2)
    if (x < 1):
        y = p_struct['min_value'] + \
                x * (p_struct['max_value'] - p_struct['min_value'])
    else:
        y = p_struct['max_value'] - \
                (x-1)*(p_struct['max_value'] - p_struct['min_value'])
        
    if (p_struct['par_mode'] == 'log'):
        y = np.power(10, y)
    
    return y