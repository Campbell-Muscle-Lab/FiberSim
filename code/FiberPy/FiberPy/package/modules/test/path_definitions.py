# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 15:50:41 2021

@author: sako231
"""

import os

# General path definitions

ROOT = os.path.dirname(__file__)
SIM_FOLDER = os.path.join(ROOT, "data")
MODEL_FILE = os.path.join(SIM_FOLDER, "sim_input/model.json")
OPTION_FILE = os.path.join(SIM_FOLDER, "sim_input/options.json")
PROTOCOL_FILE = os.path.join(SIM_FOLDER, "sim_input/pCa6_protocol.txt")
HS_STATUS_FOLDER = os.path.join(SIM_FOLDER, "sim_output/hs2")
OUTPUT_DIR = os.path.join(ROOT, "output_files")

