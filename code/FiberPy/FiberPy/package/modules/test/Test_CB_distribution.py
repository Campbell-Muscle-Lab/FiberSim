# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:00:32 2020

@author: Sarah
"""


import os,sys
import subprocess
import numpy as np
import json

import matplotlib.pyplot as plt


ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

# General path definitions

ROOT = os.path.dirname(__file__)
FIBERSIM_ROOT = os.path.join(ROOT, "..", "..","..", "..", "..", "..", "..", "FiberSim")
FIBERSIM_EXE = os.path.join(FIBERSIM_ROOT, "bin","FiberSim.exe")
INPUT_DIR = os.path.join(ROOT, "input_files")
MODEL_FILE = os.path.join(INPUT_DIR, "model_muscle_plan.json")
OPTION_FILE = os.path.join(INPUT_DIR, "options_distrib.json")
PROTOCOL_FILE = os.path.join(INPUT_DIR, "pCa9_one_time_step.txt")
#PROTOCOL_FILE = os.path.join(INPUT_DIR, "pCa45_protocol.txt")
PLOT_FOLDER = os.path.join(ROOT, "output_files")

SUBPLOTS = True


def Draw_CB_distribution(hs, filename, hsl):
    """Draw the CB distribution as a function of stretch (delta x)"""
   
    cb_extension = []
    
    cb_not_bound = []

    total_nb_cbs = len(hs.thick_fil) * hs.thick_fil[0].m_no_of_cbs

    for thick in hs.thick_fil:
        
        for cb_ind, state in enumerate(thick.cb_bound_to_a_f):
            
            if state >= 0: # CB is bound
                
                cb_x = thick.cb_x[cb_ind]
                a_f = thick.cb_bound_to_a_f[cb_ind]
                a_n = thick.cb_bound_to_a_n[cb_ind]
                
                thin_fil = hs.thin_fil[a_f]
                
                bs_x = thin_fil.bs_x[a_n]
                
                stretch = cb_x - bs_x
                #stretch = np.abs(stretch)
                
                cb_extension.append(stretch)
                
            else:
                cb_x = thick.cb_x[cb_ind]                
                a_f = thick.cb_nearest_a_f[cb_ind]
                a_n = thick.cb_nearest_a_n[cb_ind]
                
                thin_fil = hs.thin_fil[a_f]
                bs_x = thin_fil.bs_x[a_n]
                
                stretch = cb_x - bs_x
                stretch = np.abs(stretch)
                
                cb_not_bound.append(stretch)
                

                
    intervals = np.linspace(-5, 5, num=40)

    cb_distrib, inter = np.histogram(cb_extension, intervals)
    cb_distrib = cb_distrib * 100 / total_nb_cbs
    
    cb_free_distrib, inter = np.histogram(cb_not_bound, intervals)
    cb_free_distrib = cb_free_distrib / total_nb_cbs
       
    x = intervals[0:-1]
    
    if SUBPLOTS:
        
        #plt.figure()
        plt.plot(x, cb_distrib, label = f"HSL = {hsl}")
        #plt.plot(x,cb_free_distrib,'b')
        plt.ylabel("Normalized number of bound CBs (%)")
        plt.xlabel("CB stretch [nm]") 
        #plt.show() 

    
    distrib = np.asarray([x, cb_distrib])    
    np.savetxt(filename,distrib)

    
def Draw_CB_nearest_af(hs):
    """Draw the CB distribution as a function of the # nearest actin filaments"""
    
    total_nb_cbs = len(hs.thick_fil) * hs.thick_fil[-1].m_no_of_cbs
    
    cb_nearest_actin_filament = []

    for thick in hs.thick_fil:      
        for a_f in thick.cb_nearest_a_f:
            cb_nearest_actin_filament.append(a_f)
                
    intervals = range(0, len(hs.thin_fil))
    
    x_lim = len(hs.thin_fil)-1
    y_lim = 1
    
    cb_distrib, inter = np.histogram(cb_nearest_actin_filament, intervals)
    cb_distrib = cb_distrib / total_nb_cbs
        
    x = intervals[0:-1] 
    
    plt.figure()
    plt.plot(x, cb_distrib,'o')
    plt.ylabel("Normalized number of total CBs")
    plt.xlabel("# Actin filament") 
    plt.axis([0, x_lim, 0, y_lim])
    plt.show() 
    
def Draw_CB_nearest_an(hs):
    """Draw the CB distribution as a function of the # nearest actin nodes"""
    
    total_nb_cbs = len(hs.thick_fil) * hs.thick_fil[-1].m_no_of_cbs
    
    cb_nearest_actin_node = []
    
    for thick in hs.thick_fil:      
        for a_n in thick.cb_nearest_a_n:
            cb_nearest_actin_node.append(a_n)
            
    """
    cb_nearest_actin_node = []
    
    thick = hs.thick_fil[1]
    
    for a_n in thick.cb_nearest_a_n:
            cb_nearest_actin_node.append(a_n)
          
    """           
    intervals = range(0, hs.thin_fil[0].a_no_of_bs)
           
    cb_distrib, inter = np.histogram(cb_nearest_actin_node, intervals)

    cb_distrib = cb_distrib / total_nb_cbs
    
        
    x = intervals[0:-1] 

    
    distrib = np.asarray([x, cb_distrib])
    filename = os.path.join(OUTPUT_DIR, f"CB_distrib_nearest_an_{hsl}.txt")
    np.savetxt(filename,distrib)
    
    if SUBPLOTS:
        
        plt.figure()
        #plt.plot(x, cb_distrib)
        plt.bar(x,cb_distrib)
        plt.ylabel("Normalized number of total CBs")
        plt.xlabel("# Actin node") 
        #plt.xlim(300, 380)
        #plt.axis([0, x_lim, 0, y_lim])
        #plt.show() 
        
"""     

# Draw Nearest Actin Nodes distribution
    
length = [700,1000,1300]

plt.figure()

for hsl in length:
    
    with open(MODEL_FILE, "r") as jsonFile:
        data = json.load(jsonFile)

    tmp = data["muscle"]["initial_hs_length"]
    data["muscle"]["initial_hs_length"] = hsl
    
    data["thick_structure"]["m_n"] = 4
    m_n = 4

    with open(MODEL_FILE, "w") as jsonFile:
        json.dump(data, jsonFile)  
        
    with open(OPTION_FILE, "r") as jsonFile:
        data = json.load(jsonFile)

    tmp = data["options"]["log_folder"] 
    data["options"]["log_folder"] = f"output_files/dump_{hsl}"

    with open(OPTION_FILE, "w") as jsonFile:
        json.dump(data, jsonFile)  

        
    OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp_{hsl}")
    HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump_{hsl}/hs_status")
    
    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)
     
    hs_file = []
    
    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs = half_sarcomere.half_sarcomere(hs_file[-1])
 
    Draw_CB_nearest_an(hs)
    
    filename = os.path.join(OUTPUT_DIR, f"CB_distrib_nearest_an_{hsl}.txt")
    
    distribution = np.loadtxt(filename)
    node = distribution[0]
    population = distribution[1]
    
    print(population)
    

    #plt.plot(node, population, label = f"HSL = {hsl}")
    plt.bar(node, population, label = f"HSL = {hsl}")
    plt.ylabel("Normalized number of total CBs")
    plt.xlabel("# Nearest Actin Node") 
    plt.xlim(0,377)
    plt.ylim(0,0.02)
    #plt.xticks(np.arange(100, 110, 1.0))    

    
plt.legend()    
plt.title(f"m_n = {m_n}")
plt.show()
"""

# Draw filaments overlap & CB stretch distribution

length = [800]

plt.figure()

for hsl in length:
    
    with open(MODEL_FILE, "r") as jsonFile:
        data = json.load(jsonFile)

    tmp = data["muscle"]["initial_hs_length"]
    data["muscle"]["initial_hs_length"] = hsl
    
    data["thick_structure"]["m_n"] = 4
    m_n = 4

    with open(MODEL_FILE, "w") as jsonFile:
        json.dump(data, jsonFile)  
        
    with open(OPTION_FILE, "r") as jsonFile:
        data = json.load(jsonFile)

    tmp = data["options"]["log_folder"] 
    data["options"]["log_folder"] = f"output_files/dump_{hsl}"

    with open(OPTION_FILE, "w") as jsonFile:
        json.dump(data, jsonFile)  

        
    OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp_{hsl}")
    HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump_{hsl}/hs_status")
    
    cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
    subprocess.call(cmd)
     
    hs_file = []
    
    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs = half_sarcomere.half_sarcomere(hs_file[-1])
    
    # Extract the myosin and actin  x and y positions
   
    m_x = {}
    m_y = {}
    m_n = []
    
    a_x = {}
    a_y = {}
    a_n = []
    
    for thick in hs.thick_fil:
        thick_id = thick.thick_id
        m_x[thick_id] = []
        m_x[thick_id].append(thick.cb_x)
        m_y[thick_id] = []
        temp = [thick.m_y for elmt in thick.cb_x]
        m_y[thick_id].append(temp)

        
    for thin in hs.thin_fil:
        thin_id = thin.thin_id
        a_x[thin_id] = []
        a_x[thin_id].append(thin.bs_x) 
        
        a_y[thin_id] = []
        temp = [thin.a_y-0.1 for elmt in thin.bs_x]
        a_y[thin_id].append(temp)
        

    for i in range(0,len(hs.thick_fil)):
        plt.plot(m_x[i],m_y[i],'ro--')
        
    for j in range(0,len(hs.thin_fil)):
        plt.plot(a_x[j],a_y[j],'b*--')
        
    plt.xticks([])
    plt.yticks([]) 
    #plt.xlim(450,550)
    plt.title(f"HSL= {hsl} (x-y plan)")
    plt.show()
  
"""
for hsl in length:
    
    filename = os.path.join(OUTPUT_DIR, f"CB_distrib_{hsl}.txt") 
    
    HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump_{hsl}/hs_status")
    
    hs_file = []
    
    for filenames in os.listdir(HS_STATUS_FOLDER):
        filenames = os.path.join(HS_STATUS_FOLDER, filenames)

        hs_file.append(filenames)

    hs = half_sarcomere.half_sarcomere(hs_file[50])
    
    Draw_CB_distribution(hs, filename, hsl)
    plt.legend()
"""

 

     
        
   









