# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:24:28 2024

@author: Campbell
"""

import os
import json
import cv2

import numpy as np
import pandas as pd

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from natsort import natsorted

def draw_myofibril_frame(df, frame, image_file_string):
    
    fig = plt.figure( )
    