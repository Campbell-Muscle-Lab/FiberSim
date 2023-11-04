# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:44:38 2023

@author: Campbell
"""

import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from natsort import natsorted

from pathlib import Path

# Variables
top_data_folder = '../sim_data'
sim_file_string = 'sim_output/1/sim_prot_1_r1.txt'

def average_traces():
    """ Average force and Ca """
    
    # Pull off the data directories
    windows_dirs = [f for f in Path(top_data_folder).iterdir() if f.is_dir()]
    
    data_dirs = []
    for d in windows_dirs:
        data_dirs.append(str(d))
        
    data_dirs = natsorted(data_dirs)
    
    
    print(data_dirs)
    
    # for (i,d) in enumerate(data_dirs):
    #     dfs = os.path.join(os.path.abspath(str(d)),
    #                        sim_file_string)
        
    #     print(dfs)
        
    #     if (i > 10):
    #         break
        
        
if __name__ == "__main__":
    average_traces();
        
        
    
    