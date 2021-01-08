# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 09:52:23 2020

@author: Sarah
"""

import os,sys
import subprocess

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

import time

ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

# General path definitions

ROOT = os.path.dirname(__file__)
FIBERSIM_ROOT = os.path.join(ROOT, "..", "..","..", "..", "..", "..", "..", "FiberSim")
FIBERSIM_EXE = os.path.join(FIBERSIM_ROOT, "bin","FiberSim.exe")
INPUT_DIR = os.path.join(ROOT, "input_files")
MODEL_FILE = os.path.join(INPUT_DIR, "model_template.json")
OPTION_FILE = os.path.join(INPUT_DIR, "options.json")
PROTOCOL_FILE = os.path.join(INPUT_DIR, "pCa90_ramp_stretch_protocol.txt")
OUTPUT_DIR = os.path.join(ROOT, f"output_files/temp")
HS_STATUS_FOLDER = os.path.join(ROOT, f"output_files/dump/hs_status")
ERROR_FOLDER = os.path.join(ROOT, "output_files")
RESULT_FILE = os.path.join(OUTPUT_DIR, "results.txt")


t_0 = time.time()

cmd = [FIBERSIM_EXE, MODEL_FILE, OPTION_FILE, PROTOCOL_FILE, OUTPUT_DIR]
subprocess.call(cmd)

t_1 = time.time() - t_0

print(t_1)

f = plt.figure(constrained_layout=True)
f.set_size_inches([6,8])
spec = gridspec.GridSpec(nrows = 4, ncols = 1, figure = f)

ax1 = f.add_subplot(spec[0, 0])
ax1.set_ylabel('HSL (µM)')
ax1.set_xticks([])
ax2 = f.add_subplot(spec[1, 0])
ax2.set_ylabel('Force')
ax2.set_xticks([])
ax3 = f.add_subplot(spec[2, 0])
ax3.set_ylabel('Activated BS \n (proportion)')
#ax3.set_ylim([0,1])
ax4 = f.add_subplot(spec[3, 0])
ax4.set_ylabel('Attached CBs \n (proportion)')
ax4.set_xlabel('Time (s)')
#ax4.set_ylim([0,1])


fs_data = pd.read_table(RESULT_FILE)

fs_data['calcium'] = pow(10,-fs_data['pCa'])*1e+6 
       
ax1.plot('time', 'hs_length', data=fs_data, color = "tab:blue")
ax2.plot('time', 'force', data=fs_data, color = "tab:blue")
ax3.plot('time', 'a_pop_1', data=fs_data, color = "tab:blue")
ax4.plot('time', 'm_pop_1', data=fs_data, color = "tab:blue")

plt.show()

