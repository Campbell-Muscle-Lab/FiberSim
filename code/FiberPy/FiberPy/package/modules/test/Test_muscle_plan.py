# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:23:13 2020

@author: Sarah
"""

import matplotlib.pyplot as plt

import json
import subprocess

import sys, os

ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.half_sarcomere import half_sarcomere

# Choose m_n then run FiberSim

FIBERSIM_ROOT = os.path.join(ROOT, "..", "..","..", "..", "..", "..", "..", "FiberSim")
FIBERSIM_EXE = os.path.join(FIBERSIM_ROOT, "bin","FiberSim.exe")
INPUT_DIR = os.path.join(ROOT, "input_files")
MODEL_FILE = os.path.join(INPUT_DIR, "model_muscle_plan.json")
OPTION_FILE = os.path.join(INPUT_DIR, "options.json")
PROTOCOL_FILE = os.path.join(INPUT_DIR, "pCa9_one_time_step.txt")
PLOT_FOLDER = os.path.join(ROOT, "output_files")
LOG_FOLDER = os.path.join(ROOT, "output_files/dump/log_file.log")
OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp")
HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump/hs_status")

with open(MODEL_FILE, "r") as jsonFile:
    data = json.load(jsonFile)

tmp = data["muscle"]["initial_hs_length"]
data["muscle"]["initial_hs_length"] = 1615
    
data["thick_structure"]["m_n"] = 4
m_n = 4

with open(MODEL_FILE, "w") as jsonFile:
    json.dump(data, jsonFile) 

cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
subprocess.call(cmd)

# Extract the a_m_matrix from the log file

num = "hs[0]: a_m_matrix"
num_stop = "hs[0]: nearest_actin_matrix"
line_number = 0
start = 0
stop = 0

# Find the lines where the a_m_matrix is printed in the log file and save it into a 2D array

with open(LOG_FOLDER,"r") as search:            
    for line in search:
        line_number += 1
        line = line.rstrip()  # remove '/n' at end of line
        if num == line:
            start = line_number
        elif num_stop == line:
            stop = line_number-1
            break
search.close()

with open(LOG_FOLDER,"r") as f:       
    text = f.readlines() 
    
f.close()
    
temp = []
a_m_matrix = []

no_of_rows = 0
no_of_col = 0
    
for lines in text[start:stop]:
    no_of_rows += 1
    temp = [int(i) for i in lines.split()]
    no_of_col = len(temp)
    a_m_matrix.append(temp)

# Prepare the array for the thin and thick filaments (x,y) coordinates

x_a = []
y_a = []
n_a = []

x_m = []
y_m = []
n_m = []

x_mirr = []
y_mirr = []


for ind_r, row in enumerate(a_m_matrix):
    for ind_c, col in enumerate(row):
        if col < 0:
            x_a.append(ind_c)
            y_a.append(no_of_rows -1 - ind_r)
            n_a.append(col)
        elif col > 0:
            x_m.append(ind_c)
            y_m.append(no_of_rows -1 - ind_r)
            n_m.append(col)
            
        if ind_r == 0 and col != 0:
            x_mirr.append(ind_c)
            y_mirr.append(no_of_rows -1 - ind_r)
            
        if ind_c == 0 and col != 0 or ind_c == no_of_col-1 and col != 0:
            x_mirr.append(ind_c)
            y_mirr.append(no_of_rows -1 - ind_r)
          
# Plot the (x,y) map of the filaments
            
"""
            
x_a = [i + 1 for i in x_a]
y_a = [i + 1 for i in y_a]

x_m = [i + 1 for i in x_m]
y_m = [i + 1 for i in y_m]

"""
            
fig, ax = plt.subplots()
ax.scatter(x_a, y_a)
for i, txt in enumerate(n_a):
    ax.annotate(txt, (x_a[i], y_a[i]))
ax.scatter(x_m, y_m)
for i, txt in enumerate(n_m):
    ax.annotate(txt, (x_m[i], y_m[i]))

ax.scatter(x_mirr, y_mirr)

ax.set_xlabel("# column")
ax.set_ylabel("# row")
ax.set_title(f"a_m_matrix for m_n = {m_n} FIXED")

# Plot the muscle plan based on filament coordinates from the hs_status folder

hs_file = []

for filenames in os.listdir(HS_STATUS_FOLDER):
    filenames = os.path.join(HS_STATUS_FOLDER, filenames)

    hs_file.append(filenames)

hs = half_sarcomere.half_sarcomere(hs_file[-1])

m_x = []
m_y = []
m_n = []

a_x = []
a_y = []
a_n = []

for thick in hs.thick_fil:
    m_x.append(thick.m_z)
    m_y.append(thick.m_y)
    m_n.append(thick.thick_id +1)
    
for thin in hs.thin_fil:
    a_x.append(thin.a_z)
    a_y.append(thin.a_y)
    a_n.append(thin.thin_id +1)
    
    
fig, ax = plt.subplots()
ax.scatter(a_y, a_x)
for i, txt in enumerate(a_n):
    ax.annotate(txt, (a_y[i], a_x[i]))
ax.scatter(m_y,m_x)
for i, txt in enumerate(m_n):
    ax.annotate(txt, (m_y[i], m_x[i]))
    

ax.set_title("muscle plan (z and y axis inverted)")


# Check that each thin filament has a correct number of neighbouring thick filaments

actin_id = range(1, len(hs.thin_fil)+1)

for el in actin_id:
    
    counter_actin = 0
    nearest_thick = []
    
    for thick in hs.thick_fil:
        
        thick_num = thick.thick_id + 1
    
        nearest = [i + 1 for i in hs.thick_fil[thick.thick_id].nearest_actin_filaments]
    
        #print(f"Nearest actin filaments for thick filament # {thick_num} = {nearest}")
        
        if el in nearest:
            counter_actin +=1
            nearest_thick.append(thick_num)
        

    print(counter_actin)
    #print(nearest_thick)

    

    

    

    



            