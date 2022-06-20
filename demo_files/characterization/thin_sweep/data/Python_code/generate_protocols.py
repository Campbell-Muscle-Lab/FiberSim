# -*- coding: utf-8 -*-
"""
Created on Fri May  6 12:32:25 2022

@author: kscamp3
"""

import numpy as np
import pandas as pd

from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def generate_protocols():
    
    # Variables
    
    time_step = 0.001
    no_of_time_points = 700
    
    stim_start_s = 0.2
    stim_freq_Hz = 1
        
    p=dict()
    p['Ca_content'] = 1e-3
    p['k_leak'] = 6e-4
    p['k_act'] = 8.2e-2
    p['k_serca'] = 8
    p['t_open'] = 0.01
    
    protocol_file_string = '../protocols/twitch_protocol.txt'
    
    # Code
    
    # Generate a time array
    dt = time_step * np.ones(no_of_time_points)
    t = np.cumsum(dt)
    
    # Generate an activation array
    act_array = np.zeros(no_of_time_points)
    t_holder = stim_start_s
    keep_going= True
    while (keep_going):
        t_end = t_holder + p['t_open']
        if (t_end <= t[-1]):
            # There's room in the activation train for a pulse
            start_ind = int(np.round(t_holder / time_step))
            stop_ind = int(np.round(t_end / time_step))
            act_array[np.arange(start_ind,stop_ind+1)] = 1
            t_holder = t_holder + (1 / stim_freq_Hz)
        else:
            keep_going = False
    
    # Generate an array to hold the intracellular Ca
    Ca_myo = np.NaN * np.ones(no_of_time_points)
    
    def derivs(t, y):
        dy = np.zeros(y.size)
        dy[0] = (p['k_leak'] + act * p['k_act']) * y[1] - \
                    (p['k_serca'] * y[0])
        dy[1] = -dy[0]
        
        return dy
    
    y = np.array([0, p['Ca_content']])
    for i in range(no_of_time_points):
        act = act_array[i]
        sol = solve_ivp(derivs, [0, time_step], y, method='RK23')
        y = sol.y[:,-1]
        Ca_myo[i] = y[0]
        
    # Make a dataframe
    d = dict()
    d['dt'] = dt
    d['time_s'] = t
    d['Ca_myo'] = Ca_myo
    d['pCa'] = -np.log10(Ca_myo)
    d['dhsl'] = np.zeros(no_of_time_points)
    d['mode'] = -2 * np.ones(no_of_time_points)
        
    # Save Ca as a pandas array
    d = pd.DataFrame(d)
    
    # Plot
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches([4,4])
    spec=gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    ax=[]
    ax.append(fig.add_subplot(spec[0,0]))
    
    ax[0].plot(d['time_s'], d['Ca_myo'], 'b-')
    
    # Drop un-wanted columns
    d = d.drop(['time_s', 'Ca_myo'], axis=1)
    
    print(d)

    # Write to file
    d.to_csv(protocol_file_string, sep='\t', index=False)            

if __name__ == "__main__":
    generate_protocols()
