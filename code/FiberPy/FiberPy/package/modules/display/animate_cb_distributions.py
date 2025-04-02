"""@package FiberSim_analysis
Documentation for the FiberSim_analysis module

More stuff here
"""

import os
import imageio

import pandas as pd
import numpy as np

import pathlib

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from ..half_sarcomere import half_sarcomere as hs

from natsort import natsorted

def animate_cb_distributions(hs_status_directory,
                            output_gif_file,
                            frames=[]):
    """Creates an animation showing cb distributions over time """
   

    # Get hs_status_files
    hs_files = []
    
    for dirpath,_,filenames in os.walk(hs_status_directory):
        for f in filenames:
            f = os.path.abspath(os.path.join(hs_status_directory, f))
            if (f.endswith(".json")):
                hs_files.append(f)

    # Now natural sort them
    hs_files = natsorted(hs_files)

    # Limit to specific frames if required
    if (len(frames) > 0):
        temp_strings = (frames.split(':'))
        temp_indices = [int(e) for e in temp_strings]
        indices = np.arange(temp_indices[0], temp_indices[2]+1, temp_indices[1],
                            dtype=int)
        file_ends = []
        for fi in hs_files:
            p = pathlib.Path(fi)
            file_ends.append(p.parts[-1].split('_')[-1])

        hs_files_temp = []
        for ind in indices:
            test_string = '%i' % ind
            matching = [i for (i,s) in enumerate(file_ends) if test_string in s]
            if (len(matching) > 0):
                hs_files_temp.append(hs_files[matching[0]])
                
        hs_files = hs_files_temp
        
        print(hs_files)
            

    # Create an animated gif
    with imageio.get_writer(output_gif_file, mode='I') as writer:
        # Now loop through files
        for i,f in enumerate(hs_files):
            # Load the json file as an hs object
            h = hs.half_sarcomere(f)
            # Make new image file
            pre, ext = os.path.splitext(f)
            image_file_string = pre + 'cb_distrib' + '.png'
            # Draw the filament lengths
            h.draw_cb_distributions(image_file_string)
            del h
            # Add in to the gif
            image = imageio.imread(image_file_string, format='png')
            writer.append_data(image)
        

