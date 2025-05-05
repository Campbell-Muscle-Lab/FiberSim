# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:24:28 2024

@author: Campbell
"""

import os
import json
import cv2
import math
import shutil

import numpy as np
import pandas as pd

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib.patches import Rectangle

from natsort import natsorted

def draw_myofibril_frame(df, frame=1, image_file_string=[], max_t=1):
    """ Draws a myofibril """
    
    # Variables
    
    hsl_region = [0.45, 0.55]
    
    fig_hs_width = 0.5
    fig_hs_x_padding = 0.5
    fig_hs_height = 0.5
    fig_hs_y_padding = 0.5
    
    z_width = 100
    z_line_color = 'black'
    
    m_width = 50
    m_line_color = 'darkviolet'
    
    thin_y_pos = 0.9
    thin_width = 0.1
    thin_color = 'salmon'
    
    thick_y_pos = 0.5
    thick_width = 0.3
    thick_color = 'turquoise'
    
    no_of_rows = 2
    no_of_cols = 1
    
    ax_myofibril = 0
    ax_sl = 1
        
    # Work out how many half-sarcomeres there are
    command_length_cols = [x for x in d.columns.to_list() if x.endswith('command_length')]
    no_of_half_sarcomeres = len(command_length_cols)
    
    # Work out fig dimensions
    fig = plt.figure(constrained_layout=False)
    gs = fig.add_gridspec(nrows=no_of_rows, ncols=no_of_cols,
                          wspace=0.5, hspace=0.5)
    
    fig_width = (2 * fig_hs_x_padding) + (no_of_half_sarcomeres * fig_hs_width)
    fig_height = (3 * fig_hs_y_padding) + fig_hs_height
    
    fig.set_size_inches([fig_width, fig_height])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(gs[i,j]))
    
    # Make a colormap
    cm = plt.cm.viridis(np.linspace(0,1,no_of_half_sarcomeres))    
    
    # Plot half_sarcomere
    x_anchor = 0
    
    for i in range(no_of_half_sarcomeres):
        
        hsl_string = 'hs_%i_length' % (i+1)
        thin_length_string = 'hs_%i_a_length' % (i+1)
        thick_length_string = 'hs_%i_m_length' % (i+1)
        
        hs_length = d[hsl_string].iloc[frame]
        thin_length = d[thin_length_string].iloc[frame]
        thick_length = d[thick_length_string].iloc[frame]
        
        if ((i%2) == 0):
            # Z-line on left
            hs_polarity = 1
            z_line = x_anchor
            m_line = x_anchor + hs_length
        else:
            z_line = x_anchor + hs_length
            m_line = x_anchor
            hs_polarity = -1
            
        # Thin filament
        ax[ax_myofibril].add_patch(Rectangle((z_line, thin_y_pos - 0.5*thin_width),
                                              hs_polarity * thin_length,
                                                thin_width, color = thin_color))
        
        # Thick filament
        ax[ax_myofibril].add_patch(Rectangle((m_line, thick_y_pos - 0.5*thick_width),
                                             -hs_polarity * thick_length,
                                              thick_width, color = thick_color))
        
        # Z line
        ax[ax_myofibril].add_patch(Rectangle((z_line - 0.5* z_width, 0), z_width, 1,
                                             color = z_line_color))

        # M line
        ax[ax_myofibril].add_patch(Rectangle((m_line - 0.5*m_width, 0.1),
                                            m_width, 0.8, color = m_line_color))
        
        # Add label
        ax[ax_myofibril].text(0.5*(z_line + m_line), 1.1, '%i' % (i+1),
                              color=cm[i])

        # Move the anchor
        x_anchor = x_anchor + hs_length

    ax[ax_myofibril].set_xlim(0, no_of_half_sarcomeres*1500)
    
    # Now plot lengths

    t = d['time'].iloc[0:frame]
    for i in range(no_of_half_sarcomeres):
        y = d['hs_%i_length' % (i+1)][0:frame]
        ax[1].plot(t, y, color=cm[i])
        
    ax[1].set_xlim(0, max_t)
        
    
        
            
    # # Get half-sarcomere length in center
    # hs_left = math.floor(hsl_region[0] * no_of_half_sarcomeres)
    # hs_right = math.ceil(hsl_region[1] * no_of_half_sarcomeres)
    
    # # Make a new column
    # d['hsl_center'] = 0
    # for i in range(hs_left, hs_right):
    #     hsl_string = 'hs_%i_length' % i
    #     d['hsl_center'] = d['hsl_center'] + d[hsl_string]
    
        
    
    
    
    # Save
    ofs = '%s_%i.png' % (image_file_string, frame)
    fig.savefig(ofs)
    
    plt.close()
    
    
    # fig = plt.figure( )
    
def create_movie_from_image_folder(image_folder, movie_file_string):
    
    im_file_strings = [img for img in os.listdir(image_folder) if (img.endswith(".png"))]
    im_file_strings = natsorted(im_file_strings)
    
    print(im_file_strings)
    
    frame = cv2.imread(os.path.join(image_folder, im_file_strings[0]))
    height, width, layers = frame.shape
    
    video = cv2.VideoWriter(movie_file_string, -1, 30, (width, height))
    
    for im in im_file_strings:
        video.write(cv2.imread(os.path.join(image_folder, im)))
        
    # cv2.destroyAllWindows()
    video.release()
    
    
if __name__ == "__main__":
    fs = 'd:/ken/github/campbellmusclelab/models/fibersim/demo_files/myofibrils/spocs_3/sim_data/sim_output/5/sim_pCa_57_s_0_r1.txt'
    
    d = pd.read_csv(fs, sep='\t')
    
    od = 'd:/temp/test_images/ken'
    od2 = 'd:/temp/test_images'
    mfs = 'd:/temp/myofibril.mp4'
    
    print(d.shape)
    
    max_t = d['time'].iloc[-1]
    
    # Clean the image directory
    try:
        print('Trying to remove %s' % od2)
        shutil.rmtree(od2, ignore_errors=True)
    except OSError as e:
        print('Error: %s: %s' % (od2, e.strerror))
        
    if not os.path.isdir(od2):
        os.makedirs(od2)
    
    for i in np.arange(1, d.shape[0], 100):
        draw_myofibril_frame(d, i, od, max_t)
    
    create_movie_from_image_folder(od2, mfs)