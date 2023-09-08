# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 12:35:35 2023

@author: Campbell
"""

import os
import sys
import json

from pathlib import Path

def analyze_status_file(status_file = ''):
    
    status_file = '../sim_data/force_velocity/isotonic/sim_output/1/status_dumps_1_r1/hs_1_time_step_480.json'
    
    with open(status_file, 'r') as f:
        hs = json.load(f)
        
    print(hs)
    
if __name__ == "__main__":
    analyze_status_file()
    
    
    