# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 13:23:35 2021

@author: kscamp3
"""

import json

import pandas as pd
import numpy as np

def create_frames_file():
    """ Generates a frames file """

    frames_file_string = '../frames/frames2.json'
    result_file_string = '../results/results.txt'

    # Load up FiberSim results
    d = pd.read_csv(result_file_string, sep='\t')
    
    hsl = d['hs_length'].to_numpy()
    hsl_ref = hsl[8000]
    print(hsl_ref)
    offset = hsl_ref - 609

    frame = []

    for i in range(1,10000,20):
        
        fr = dict()
        fr['status_file'] = 'snapshots/hs_1_time_step_%i.json' % i
        fr['image_file'] = 'renders/hs_c_zone_%i.png' % i
        c = dict()
        loc = dict()
        loc['x'] = hsl[i] - offset
        loc['y'] = 83
        loc['z'] = 33.4
        c['location'] = loc
        rot = dict()
        rot['x'] = 81.1
        rot['y'] = 0.8
        rot['z'] = -107
        c['rotation'] = rot
        fr['camera'] = c

        frame.append(fr)

    for i in range(1,10000,20):

        fr = dict()
        fr['status_file'] = 'snapshots/hs_1_time_step_%i.json' % i
        fr['image_file'] = 'renders/hs_i_band_%i.png' % i
        c = dict()
        loc = dict()
        loc['x'] = 94
        loc['y'] = 43.3
        loc['z'] = 33.2
        c['location'] = loc
        rot = dict()
        rot['x'] = 85.9
        rot['y'] = 0.836
        rot['z'] = -94.7
        c['rotation'] = rot
        fr['camera'] = c

        frame.append(fr)

    frames = dict()
    frames['frames'] = frame
    
    with open(frames_file_string, 'w') as f:
        json.dump(frames, f, indent=4)

    print('no_of_frames: %i' % len(frames['frames']))

if __name__ == "__main__":
    create_frames_file()
    
