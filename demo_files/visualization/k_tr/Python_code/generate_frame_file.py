# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 15:00:19 2021

@author: kscamp3
"""

import json

import numpy as np

def generate_frame_file():
    # Generates a list of frames

    template_frame_file_string = '../frames/frames_template.json'
    base_status_file = 'sim_output/hs/hs_1_time_step'
    base_render_file = 'renders/im'
    frames_file = '../frames/frames_ktr.json'
    
    with open(template_frame_file_string, 'r') as in_f:
        templ = json.load(in_f)
    print(templ)

    n_steps = np.asarray([75, 750, 75, 75])
    
    out = []
    no_of_anchors = len(templ['frames'])
    
    file_indices = np.arange(1, 10000, 10)
    fc = 0
    
    for i in range(no_of_anchors-1):
        start = templ['frames'][i]
        stop = templ['frames'][i+1]

        # Create a camera struct
        cam = []

        # Interpolate the positions and rotations
        for j in range(n_steps[i]):
            camera = dict()
            
            for t in ['location', 'rotation']:
                camera[t] = dict()
                for c in ['x', 'y', 'z']:
                    a = start['camera'][t][c]
                    b = stop['camera'][t][c]
                    p = np.linspace(a,b, n_steps[i])
                
                    camera[t][c] = p[j]
            
            sf = ('%s_%.0f.json' % (base_status_file, file_indices[fc]))
            rf = ('%s_%.0f.png' % (base_render_file, file_indices[fc]))
            out_f = dict()
            out_f['status_file'] = sf
            out_f['image_file'] = rf
            out_f['camera'] = camera
            
            out.append(out_f)
            
            fc = fc+1
        
        out_frames = dict()
        out_frames['frames'] = out
            
        with open(frames_file, 'w') as outfile:
            json.dump(out_frames, outfile, indent=4)
                    
        #                       x1 = start['camera']['location']['x']
        # y1 = start['camera']['location']['y']
        # z1 = start['camera']['location']['z']
        
        # cx1 = start['camera']['rotation']['x']
        # cy1 = start['camera']['rotation']['y']
        # cz1 = start['camera']['rotation']['z']
        
        
    print(fc)

if __name__ == "__main__":
    generate_frame_file()