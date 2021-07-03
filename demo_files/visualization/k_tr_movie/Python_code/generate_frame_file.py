# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 15:00:19 2021

@author: kscamp3
"""

import json

import numpy as np

def generate_frame_file():
    # Generates a list of frames

    frame_indices = np.arange(3500, 8000, 10)
    
    cam_loc_x = 875
    cam_loc_y = 30
    cam_loc_z = 17.32
    
    cam_rot_x = 90
    cam_rot_y = 0
    cam_rot_z = 90

    base_status_file = 'sim_output/hs/hs_1_time_step'
    base_render_file = 'renders/im'
    frames_file = '../frames/frames_ktr.json'

    out=[]

    for ind in frame_indices:

        camera = dict()
        
        for t in ['location', 'rotation']:
            camera[t] = dict()
        
        camera['location']['x'] = cam_loc_x
        camera['location']['y'] = cam_loc_y
        camera['location']['z'] = cam_loc_z
        
        camera['rotation']['x'] = cam_rot_x
        camera['rotation']['y'] = cam_rot_y
        camera['rotation']['z'] = cam_rot_z
        
        sf = ('%s_%.0f.json' % (base_status_file, ind))
        rf = ('%s_%.0f.png' % (base_render_file, ind))
        out_f = dict()
        out_f['status_file'] = sf
        out_f['image_file'] = rf
        out_f['camera'] = camera
        
        out.append(out_f)

    out_frames = dict()
    out_frames['frames'] = out
        
    with open(frames_file, 'w') as outfile:
        json.dump(out_frames, outfile, indent=4)

if __name__ == "__main__":
    generate_frame_file()