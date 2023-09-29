# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 18:38:37 2022

@author: kscamp3
"""

import os

import numpy as np
import pandas as pd

from scipy.integrate import solve_ivp

def write_protocol_to_file(prot, prot_file_string):
    """ Writes a protocol defined as a Pandas dataframe to a file_string
        defined with an absolute path """
    
    # Parse the file_string, and check the parent folder exists
    # If not, create the directory
    parent_dir = os.path.dirname(prot_file_string)
    
    if not os.path.isdir(parent_dir):
        print('Creating parent dir: %s' % parent_dir)
        os.makedirs(parent_dir)
    
    # Write prot to file
    print('Writing protocol to: %s' % prot_file_string)
    prot.to_csv(prot_file_string, index=False, sep='\t',
                float_format="%.5f")

def create_length_control_protocol(time_step=0.0001,
                                   n_points=1000,
                                   initial_pCa= 9.0,
                                   step_pCa = 4.5,
                                   step_up_pCa_s = 0.10,
                                   step_down_pCa_s = [],
                                   mode_vector = [],
                                   delta_hsl=[]):
    """ Creates a length control protocol as a Pandas dataframe """
    
    dt = time_step * np.ones(n_points)
    t = np.cumsum(dt)
    pCa = initial_pCa * np.ones(n_points)
    pCa[t > step_up_pCa_s] = step_pCa
    if (not (step_down_pCa_s == [])):
        pCa[ t > step_down_pCa_s] = initial_pCa
    
    # Set dhsl as 0 if undefined (thus isometric)
    # Otherwise, fill based on delta_hsl
    if (delta_hsl == []):
        dhsl = np.zeros(n_points)
    else:
        dhsl = delta_hsl
    
    # Set mode_vector as -2 if undefined (thus length control)
    if (mode_vector == []):
        mode = -2 * np.ones(n_points)
    else:
        mode = mode_vector
           
    
    # Assemble
    d = {'dt': dt, 'pCa': pCa, 'dhsl': dhsl, 'mode': mode}
    df = pd.DataFrame(data=d)
    
    return df

def create_force_control_protocol(time_step=0.0001, n_points=1000,
                              initial_pCa=9.0, step_pCa=4.5,
                              step_pCa_s=0.01,
                              mode_vector=[],
                              iso_f=0,
                              iso_start_s=0.08):
    """ Creates a force control protocol as a Pandas dataframe """
    
    dt = time_step * np.ones(n_points)
    t = np.cumsum(dt)
    pCa = initial_pCa * np.ones(n_points)
    pCa[t > step_pCa_s] = step_pCa
    dhsl = np.zeros(n_points)
    if (mode_vector == []):
        mode = -2 * np.ones(n_points)
        mode[t > iso_start_s] = iso_f
    else:
        mode = mode_vector
    
    # Assemble
    d = {'dt': dt, 'pCa': pCa, 'dhsl': dhsl, 'mode': mode}
    df = pd.DataFrame(data=d)
    
    return df

def create_twitch_protocol(time_step=0.001, n_points=600,
                           stimulus_times_s = [0.1],
                           Ca_content=1e-3, stimulus_duration_s=0.01,
                           k_leak=6e-4, k_act=8.2e-2, k_serca=20,
                           dhsl = [],
                           mode_vector=[]):
    """ Creates a twitch protocol """
    
    dt = time_step * np.ones(n_points)
    activation = np.zeros(n_points)
    myofil_Ca = np.nan * np.ones(n_points)
    
    for st in stimulus_times_s:
        act_start_index = round(st / time_step)
        act_stop_index = round((st + stimulus_duration_s) / time_step)
        activation[act_start_index:(act_stop_index+1)] = 1
    
    # Solve two compartment problem
    
    def derivs(t, y):
        dydt = np.zeros(2)
        dydt[0] = ((k_leak + (act * k_act)) * y[1]) - (k_serca * y[0]);
        dydt[1] = -dydt[0]
        return dydt
    
    y_temp = [0, Ca_content]
    for i in range(n_points):
        act = activation[i]
        sol = solve_ivp(derivs, [0, time_step], y_temp, method='RK23')
        y_temp = sol.y[:,-1]
        myofil_Ca[i] = y_temp[0]
    
    pCa = -np.log10(myofil_Ca)
    
    if (mode_vector == []):
        mode_vector = -2 * np.ones(n_points)
    
    if (dhsl == []):
        dhsl = np.zeros(n_points)
    
    # Assemble
    d = {'dt': dt, 'pCa': pCa, 'dhsl': dhsl, 'mode': mode_vector}
    df = pd.DataFrame(data=d)

    return df    
    
    
    
    
    


