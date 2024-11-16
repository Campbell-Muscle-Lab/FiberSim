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

from matplotlib.patches import Rectangle

from natsort import natsorted

def draw_myofibril_frame(df, frame=1, image_file_string=[]):
    """ Draws a myofibril """
    
    # Variables
    fig_hs_width = 0.5
    fig_hs_x_padding = 0.5
    fig_hs_height = 0.5
    fig_hs_y_padding = 0.5
    
    ax_myofibril = 0
    
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
    
        
    # Work out how many half-sarcomeres there are
    command_length_cols = [x for x in d.columns.to_list() if x.endswith('command_length')]
    no_of_half_sarcomeres = len(command_length_cols)
    
    # Work out fig dimensions
    fig = plt.figure(constrained_layout=False)
    gs = fig.add_gridspec(nrows=1, ncols=1)
    
    fig_width = (2 * fig_hs_x_padding) + (no_of_half_sarcomeres * fig_hs_width)
    fig_height = (2 * fig_hs_y_padding) + fig_hs_height
    
    fig.set_size_inches([fig_width, fig_height])
    ax = []
    ax.append(fig.add_subplot(gs[0,0]))
    
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
        else:
            x_anchor = x_anchor + hs_length
            hs_polarity = -1

        # Z line
        ax[ax_myofibril].add_patch(Rectangle((x_anchor - 0.5* z_width, 0), z_width, 1,
                                             color = z_line_color))

        # M line
        ax[ax_myofibril].add_patch(Rectangle((x_anchor + hs_polarity*hs_length - 0.5*m_width, 0.1),
                                            m_width, 0.8, color = m_line_color))
        
        # Thin filament
        ax[ax_myofibril].add_patch(Rectangle((x_anchor, thin_y_pos - 0.5*thin_width),
                                              hs_polarity * thin_length,
                                               thin_width, color = thin_color))
        
        # Thick filament
        ax[ax_myofibril].add_patch(Rectangle((x_anchor + hs_polarity*hs_length, thick_y_pos - 0.5*thick_width),
                                             hs_polarity * thick_length,
                                             thick_width, color = thick_color))

        # Move the anchor
        if (hs_polarity == 1):        
            x_anchor = x_anchor + hs_length
        
    plt.xlim(0, 70000)
    plt.ylim(-1, 2)
    
    
    
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
    fs = 'd:/ken/github/campbellmusclelab/publications/paper_fibersim_sequential_relaxation/simulations/kens_playground/f/sim_data/sim_output/1/sim_pCa_59_s_0_r1.txt'
    
    d = pd.read_csv(fs, sep='\t')
    
    od = 'd:/temp/test_images/ken'
    od2 = 'd:/temp/test_images'
    mfs = 'd:/temp/myofibril.mp4'
    
    for i in np.arange(1, 4500, 10):
        draw_myofibril_frame(d, i, od)
    
    create_movie_from_image_folder(od2, mfs)