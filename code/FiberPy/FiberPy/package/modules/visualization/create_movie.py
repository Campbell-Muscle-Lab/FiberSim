# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 15:46:10 2021

@author: sako231
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': "Arial"})

import pandas as pd
from natsort import natsorted

import cv2

from package.modules.half_sarcomere import half_sarcomere

import matplotlib.gridspec as gridspec

def create_hs_movie(dump_folder, movie_folder):
    """
    Creates an animated movie showing pCa, force and the filament states
    """
    
    # extracting the dump files
    
    hs_file = []

    for filenames in os.listdir(dump_folder):
        filenames = os.path.join(dump_folder, filenames)

        hs_file.append(filenames)
        
    # naturally sorting the dump files

    hs_file = natsorted(hs_file)
    
    # create the force and pCa arrays
    
    force_calcium_data = get_hs_data(hs_file)
    
    # create the snapshots array
    
    frames = []
    
    for i, file in enumerate(hs_file):
        
        print('Rendering snapshot of "%s"' % file) 
        
        hs = half_sarcomere.half_sarcomere(file)
        
        if not os.path.isdir(movie_folder):
            os.makedirs(movie_folder)
        
        snapshot_folder = os.path.join(movie_folder, "snapshots")
        
        if not os.path.isdir(snapshot_folder):
            os.makedirs(snapshot_folder)
        
        snapshot = os.path.join(snapshot_folder, f"hs_snapshot_{i}.png")
        
        frames.append(create_hs_snapshot(hs, force_calcium_data, i, snapshot)) # create array of image files
        
    # create the movie
    
    img_array = []
    for f in frames:
        im = cv2.imread(f)
        height, width, layers = im.shape
        size = (width,height)
        img_array.append(im)
        
    video_path = os.path.join(movie_folder, "snapshots_movie.mp4")
    
    out = cv2.VideoWriter(video_path, cv2.VideoWriter_fourcc(*'mp4v'), 10, size)
    
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()
    print('Movie created')
        
        
def create_hs_snapshot(hs, pd_data, time_step, im_file_str):

    # Get filaments data
    
    ### Thick Filaments ###
    
    for i, t in enumerate(hs['thick']):
        if (i==0):
            
            # CB states
            m_1 = np.zeros([len(hs['thick']),len(t['cb_state'])])
            m_2 = np.zeros([len(hs['thick']),len(t['cb_state'])])
            m_3 = np.zeros([len(hs['thick']),len(t['cb_state'])])
            
            # CB x-positions
            
            m_x = np.zeros([len(hs['thick']),len(t['cb_state'])])
            
            # Thick filament
            
            thick = np.zeros(len(t['cb_state']))
            
        m_x[i,:] = t['cb_x'] 
        
        m_1[i,:] = [(x == 1) for x in t['cb_state']] # SRX
        m_2[i,:] = [(x == 2) for x in t['cb_state']] # DRX
        m_3[i,:] = [(x == 3) for x in t['cb_state']] # FG
        
    # Averaging over thick filaments
    
    m_1 = np.mean(m_1, 0)
    m_2 = np.mean(m_2, 0)
    m_3 = np.mean(m_3, 0)
    
    m_x = np.mean(m_x, 0)
    
    # Averaging for each node
    
    size = len(hs["thick"][0]["node_forces"])
    m_1_a = np.zeros(size)
    m_2_a = np.zeros(size)
    m_3_a = np.zeros(size)
    m_x_a = np.zeros(size)
  
    for i, node in enumerate(hs["thick"][0]["node_forces"]):
        
        idx = i * 6
        
        m_1_a[i] = np.mean(m_1[idx:idx+6])
        m_2_a[i] = np.mean(m_2[idx:idx+6])
        m_3_a[i] = np.mean(m_3[idx:idx+6])
        
        m_x_a[i] = np.mean(m_x[idx:idx+6])

    
    # Thick filament -> add bare-zone
            
    bare_x = m_x_a[0] + 13.5 + hs['thick'][0]["m_lambda"] # cb_x_0 + m_rl + lambda 
    
    thick_x = np.insert(m_x, 0, bare_x)
    thick = np.insert(thick, 0 , 0)
    
            
    ### Thin filaments ###
    
    for i, t in enumerate(hs['thin']):
        if (i==0):
            
            # BS states
            bs_1 = np.zeros([len(hs['thin']),len(t['bs_state'])])
            bs_2 = np.zeros([len(hs['thin']),len(t['bs_state'])])
            
            # BS x-positions
            
            bs_x = np.zeros([len(hs['thin']),len(t['bs_state'])])
            
            # Thin filament
            
            thin = np.zeros(len(t['bs_state'])) + 2
            
        bs_x[i,:] = t['bs_x'] 
        
        bs_1[i,:] = [(x == 1) for x in t['bs_state']]   # OFF
        bs_2[i,:] = [(x == 2) for x in t['bs_state']] # ON

    # Averaging over thin filaments
    
    bs_1 = np.mean(bs_1, 0)
    bs_2 = np.mean(bs_2, 0)
    
    bs_x = np.mean(bs_x, 0)
    
    # Averaging for each node
    
    size = int(len(hs["thin"][0]["bs_state"])/2)
    
    bs_1_a = np.zeros(size)
    bs_2_a = np.zeros(size)

    bs_x_a = np.zeros(size)
    
    
    for i in range(0, size):
        
        idx = i * 2
        
        bs_1_a[i] = np.mean(bs_1[idx:idx+2])
        bs_2_a[i] = np.mean(bs_2[idx:idx+2])
        
        bs_x_a[i] = np.mean(bs_x[idx:idx+2])
        
        
    ### Finding the D, C, P and M zones

    
    pc_node_idx_0 = hs['thick'][0]['pc_node_index'][0]
    
    pc_idx_0 = pc_node_idx_0 * hs['thick'][0]['m_cbs_per_node']
    
    pc_x_0 = m_x[pc_idx_0] # average x-position of the first C-protein
    
    pc_node_idx_1 = hs['thick'][0]['pc_node_index'][-1]
    
    pc_idx_1 = pc_node_idx_1 * hs['thick'][0]['m_cbs_per_node']
    
    pc_x_1 = m_x[pc_idx_1] # average x-position of the last C-protein
    
    ### PLOTS ###

    
    # fig, ax = plt.subplots(5, figsize=(8, 6), dpi=200)
    # fig.subplots_adjust(hspace=0.7)

    # fig = plt.figure(constrained_layout=True, figsize=(12, 8), dpi = 200)
    # spec = fig.add_gridspec(11, 1, hspace=0.05, height_ratios = [1,1,1,1,1,1,1,2,2,1,1])
    
    # spec.update(hspace=0.05)
    
    fig = plt.figure(constrained_layout=True, figsize=(3.5, 4), dpi = 200)
    spec = fig.add_gridspec(6, 1, hspace=0, wspace =0, height_ratios = [1,1,0.8,1.6,2,2.2])
    
    spec.update(hspace=0.05)

    
    # for i in range(5):
    #     ax[i] = fig.add_subplot(spec[i,0])
    
    ax = []
    
    ax0 = fig.add_subplot(spec[0,0])
    ax1 = fig.add_subplot(spec[1,0])
    ax2 = fig.add_subplot(spec[3,0])
    ax3 = fig.add_subplot(spec[4,0])
    ax4 = fig.add_subplot(spec[5,0])
    
    ax = [ax0, ax1, ax2, ax3, ax4]
        
    for i in range(5):
        
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['bottom'].set_visible(False)
        ax[i].set_xticks([])    
        ax[i].set_xlim([-5, hs["hs_data"]["hs_length"] + 10])
        
        for axis in ['top','bottom','left','right']:
            ax[i].spines[axis].set_linewidth(1.5)
            
        for tick in ax[i].yaxis.get_major_ticks():
            tick.label.set_fontsize(10) 
            
        ax[i].tick_params(direction = "out", length = 6, width = 1.5)

    ax[4].spines['bottom'].set_visible(True)
    
    # 0) pCa plot
    
    ax[0].plot('time', 'pCa', data = pd_data[0:time_step+1], color = "black", linewidth = 2)  
    
    # Y axis
    
    ax[0].set_ylim([4,9.5])
    ax[0].spines['left'].set_bounds(4,9)
    ax[0].set_yticks([4,9]) 
    ax[0].text(-0.2, 7, "pCa" , fontsize = 10, transform = ax[0].transData, ha='center', va='center')
    ax[0].invert_yaxis()
    
    # X axis
    
    #max_time = max(pd_data["time"]) 
    
    max_time = multiple_greater_than(max(pd_data["time"]),
               0.1*np.power(10, np.ceil(np.log10(max(pd_data["time"])))))
    
    ax[0].set_xlim(0, max_time)

    # 1) Force plot
    
    ax[1].plot('time', 'force', data = pd_data[0:time_step+1], color = "black", linewidth = 2) 

    # Y axis
    
    max_y_axis = multiple_greater_than(max(pd_data["force"]),
               0.1*np.power(10, np.ceil(np.log10(max(pd_data["force"])))))
       
    ax[1].text(-0.2, 0.4*max_y_axis, "Force\n(kN $\mathregular{m^{-2}}$)" , fontsize = 10, 
               ha='center', va='center', transform = ax[1].transData, linespacing=0.9)
    ax[1].spines['bottom'].set_visible(True)
    ax[1].spines['bottom'].set_bounds(0, max_time)
    ax[1].spines['left'].set_bounds(0, max_y_axis)
    
    
    y_ticks = [0, max_y_axis]
    ax[1].set_ylim(y_ticks)
    ax[1].set_yticks(y_ticks)
    

    ax[1].set_xlim([-0.005, max(pd_data["time"])]) 
    ax[1].set_xticks([0, max_time]) 
    ax[1].set_xlabel("Time (s)", fontsize = 10, labelpad = -10)
    
    for tick in ax[1].xaxis.get_major_ticks():
        tick.label.set_fontsize(10) 
    
    # X axis 
    
    #ax[1].set_xlim(0, max(pd_data["time"]))
    
    # 2) Thin filament states
    
    ax[2].plot(bs_x_a, bs_1_a, '-', linewidth = 2, color = "tab:red")
    ax[2].plot(bs_x_a, bs_2_a, '-', linewidth = 2, color = "tab:green")
    
    # Y axis 
    
    ax[2].set_ylim([-0.1,1.1])
    ax[2].spines['left'].set_bounds(0, 1)
    ax[2].set_yticks([0, 1])
    
    ax[2].text(-360, 0.5, "Thin \n relative \n population" , fontsize = 10, 
               ha='center', va='center', transform = ax[2].transData, linespacing=0.9)
    
    
    # Label curves
    
    ax[2].text(bs_x_a[-1] + 20, bs_1_a[-1] - 0.1, "Inactive", fontsize = 10)
    
    ax[2].text(bs_x_a[-1] + 20, bs_2_a[-1] - 0.1, "Active", fontsize = 10)

    
    # 3) Thin and thick filaments
    
    ax[3].spines['left'].set_visible(False)
    ax[3].set_yticks([])
    ax[3].set_ylim([-0.5,3])
    
    ax[3].plot(bs_x, thin, color="tab:grey", linewidth = 6, solid_capstyle="butt")
    ax[3].plot(thick_x, thick, color="tab:grey", linewidth = 14, alpha = 0.5, solid_capstyle="butt")
    
    # Labeling the Z and M line
    
    ax[3].axvline(x=0,ymin=0,ymax=1, linestyle = '-', linewidth = 1, color = 'black')
    ax[3].text(-250, 1, "Z-line", fontsize = 10)
    
    ax[3].axvline(x=hs["hs_data"]["hs_length"],ymin=0,ymax=1, linestyle = '-', linewidth = 4, color = 'black', alpha = 0.7)
    ax[3].text(hs["hs_data"]["hs_length"] + 40, 1, "M-line", fontsize = 10)
    
    # Labeling the different sarcomere zones        
    
    x_D = m_x[-1] + (abs(m_x[-1] - pc_x_1))/2
    
    ax[3].text(x_D - 15, -0.3, "D", fontsize = 9)
    
    x_C = pc_x_1 + (abs(pc_x_0 - pc_x_1))/2
    
    ax[3].text(x_C - 15, -0.3, "C", fontsize = 9)
    
    x_P = pc_x_0 + (abs(m_x[0] - pc_x_0))/2
    
    ax[3].text(x_P - 15 , -0.3, "P", fontsize = 9)
    
    x_M = m_x[0] + (abs(m_x[0] - bare_x))/2
    
    ax[3].text(x_M - 25, -0.3, "B", fontsize = 9)
    
    ax[3].axvline(x=m_x[-1],ymin=0,ymax=0.3, linestyle = '-', linewidth = 2, color = 'black', solid_capstyle="butt")
    
    ax[3].axvline(x=pc_x_1,ymin=0,ymax=0.3, linestyle = '-', linewidth = 2, color = 'black', solid_capstyle="butt")
    
    ax[3].axvline(x=pc_x_0,ymin=0,ymax=0.3, linestyle = '-', linewidth = 2, color = 'black', solid_capstyle="butt")
    
    ax[3].axvline(x=m_x[0],ymin=0,ymax=0.3, linestyle = '-', linewidth = 2, color = 'black', solid_capstyle="butt")
    
    # Labeling the filaments
    
    ax[3].axvline(x=0,ymin=0,ymax=1, linestyle = '-', linewidth = 6, color = 'black')
    ax[3].text(400 , 0.65, "Thick filament", fontsize = 10)
    
    ax[3].axvline(x=0,ymin=0,ymax=1, linestyle = '-', linewidth = 6, color = 'black')
    ax[3].text(50, 2.35, "Thin filament", fontsize = 10)
    
    # 4) Thick filament states
    
            
    ax[4].plot(m_x_a, m_1_a, '-', linewidth = 2, color = "#4d4d4aff", label = "SRX")
    ax[4].plot(m_x_a, m_2_a, '-', linewidth = 2, color = "#00297ffd", label = "DRX")
    ax[4].plot(m_x_a, m_3_a, '-', linewidth = 2, color = "#a0001bff", label = "FG")
    
    # Label curves
    
    # ax[4].text(m_x_a[0] + 10, m_1_a[0] - 0.1, "SRX", fontsize = 10)
    
    # ax[4].text(m_x_a[0] + 10, m_2_a[0] - 0.1, "DRX", fontsize = 10)
    
    # ax[4].text(pc_x_0 + 40, m_3_a[pc_node_idx_0] + 0.05, "FG", fontsize = 10)
    
    ax[4].legend(loc = [0.95,0.2], frameon = False, fontsize = 8)
    
    # Y axis
    
    ax[4].set_ylim([-0.5,1.1])
    ax[4].spines['left'].set_bounds(0, 1)
    ax[4].set_yticks([0, 1])
    ax[4].spines['left'].set_position(("data", m_x_a[-1] - 5))
    
    ax[4].text(145, 0.55, "Thick \n relative \n population" , fontsize = 10, 
               ha='center', va='center', transform = ax[4].transData, linespacing=0.9)
    
    # X axis 
    
    ax[4].spines['bottom'].set_bounds(0, hs["hs_data"]["hs_length"])
    
    ax[4].set_xlim([-5, hs["hs_data"]["hs_length"] + 10]) 
    ax[4].set_xticks([0, hs["hs_data"]["hs_length"]]) 
    ax[4].set_xlabel("Distance from Z-line \n (nm)", fontsize = 10, labelpad = -10)
    
    for axis in ['top','bottom','left','right']:
        ax[4].spines[axis].set_linewidth(1.5)
            
    for tick in ax[4].yaxis.get_major_ticks():
        tick.label.set_fontsize(10) 
        
    for tick in ax[4].xaxis.get_major_ticks():
        tick.label.set_fontsize(10) 
    

    # Image output
    # print('Saving snapshot to "%s"' % im_file_str) 
    
    fig.savefig(im_file_str, dpi=200, bbox_inches = "tight")
    plt.close()
    
    return im_file_str

    
def get_hs_data(dump_files):
    
    """
    Extract force and pCa for each time-step, and store into dataframe
    """
    
    force_data = []
    pCa_data = []
    time_data = []

    for file in dump_files:
        
        hs = half_sarcomere.half_sarcomere(file)
        
        force_data.append(hs["hs_data"]["hs_force"])
        pCa_data.append(hs["hs_data"]["pCa"])
        time_data.append(hs["hs_data"]["time"])
        
    # create dataframe
    
    d = {'time': time_data, 'pCa': pCa_data, 'force' : force_data}
    
    df = pd.DataFrame(data = d)
    
    df["force"] /= 1000
    
    return df

def multiple_greater_than(val, mult):
    # Returns a multiple greater than

    return (mult * np.ceil(val/mult))
