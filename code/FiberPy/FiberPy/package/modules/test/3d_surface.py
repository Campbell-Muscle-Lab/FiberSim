# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:59:05 2021

@author: sako231
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:08:32 2021

@author: sako231
"""

import os
import json

import path_definitions as pd
import kinetics_data as kd

import numpy as np

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm
import math


import pandas

def compute_m_kinetics_3D_plots():
     
    # load calculated rates

    outfile = os.path.join(pd.OUTPUT_DIR, "calc_dict.json")  
    with open(outfile, 'r') as f:
        calculated_rates_dict = json.load(f)
        
    ## 3D plots
    
    angle_array = np.arange(kd.A_MIN,kd.A_MAX, kd.A_STEP)

    for i, iso in enumerate(calculated_rates_dict):
        
        rates_file = os.path.join(pd.OUTPUT_DIR, f"rate_equations_iso_{i}.csv")  
        rate_data = pandas.read_csv(rates_file)
        
        for j, trans in enumerate(iso):
            
            if trans["trans_type"] == 'a':
                
                idx = trans["trans_idx"]
                rate_type = trans["rate_type"]

                calc_rate = np.transpose(np.array(trans["calculated_rate"]))
                
                stretches = np.array(rate_data["stretch"] + kd.X_STEP/2)
                angles = angle_array
                
                X = stretches
                Y = angles + kd.A_STEP/2
                X, Y = np.meshgrid(X, Y)
                Z = calc_rate
                
                real_rate = Z - Z
                for af in range(0, len(angle_array)):
                    
                    align_factor = abs(-math.cos((angle_array[af] + kd.A_STEP/2)*math.pi/180))
                    real_rate[af, :] = np.array(rate_data[f"Transition # {idx} ({rate_type})"]) * align_factor #transposed!
        
                
                # set up a figure twice as wide as it is tall
                fig = plt.figure(figsize=plt.figaspect(0.5), dpi = 300)
                ax = fig.add_subplot(1, 2, 1, projection='3d')
                #ax.set_zlim(0,100)
                ax.set_xlim(-8,8)
                ax.set_ylim(80,190)
                #ax.contour3D(X, Y, Z, 10, zdir='y', rstride=1, cstride=1, cmap='plasma')
                surf = ax.plot_surface(X, Y , real_rate - Z, rstride=1, cstride=1, color = "tab:grey", alpha=1)
                #ax.scatter3D(X, Y, Z, '-o')

                ax.view_init(20,135)
                
                fig = plt.figure(figsize=plt.figaspect(0.5), dpi = 300)
                ax = fig.add_subplot(1, 2, 1, projection='3d')
                ax.set_zlim(0,200)
                ax.set_xlim(-8,8)
                ax.set_ylim(80,190)
                #ax.contour3D(X, Y, Z, 10, zdir='y', rstride=1, cstride=1, cmap='plasma')
                ax.plot_surface(X, Y , Z, rstride=1, cstride=1, color = "blue", alpha=0.5)
                ax.plot_surface(X, Y , real_rate, rstride=1, cstride=1, color = "black", alpha=0.5)
                #ax.scatter3D(X, Y, Z, '-o')

                ax.view_init(0,240) # (0,60) or (0, 240)
                
                fig = plt.figure(figsize=plt.figaspect(0.5), dpi = 300)
                ax = fig.add_subplot(1, 2, 1, projection='3d')
                ax.set_zlim(0,200)
                ax.set_xlim(-8,8)
                ax.set_ylim(80,190)
                #ax.contour3D(X, Y, Z, 10, zdir='y', rstride=1, cstride=1, cmap='plasma')
                ax.plot_surface(X, Y , Z, rstride=1, cstride=1, color = "blue", alpha=0.5)
                ax.plot_surface(X, Y , real_rate, rstride=1, cstride=1, color = "black", alpha=0.5)
                #ax.scatter3D(X, Y, Z, '-o')

                ax.view_init(0,60) # (0,60) or (0, 240) or (0, 120) or (0, 180)/(0,0)

            
compute_m_kinetics_3D_plots()
