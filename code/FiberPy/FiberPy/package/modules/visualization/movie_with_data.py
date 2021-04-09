# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 22:31:21 2021

@author: ken
"""

import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class movie_with_data():
    """ Class to create a movie that includes screen-shots and plots """

    def __init__(self, blah):
        """ Constructor """
        
        self.initialise_figure()
        
        fs = 'C:/ken/GitHub/CampbellMuscleLab/models/FiberSim/demo_files/visualization/movie_with_data/renders/hs_10.png'
        self.display_snapshot(fs)
        
        self.save_fig_as_image('c:/temp/ken.png')
    
    def initialise_figure(self):
        """ Initialises figure """

        # Make a figure iwth 1 large panel to show images and
        # a column of data panels on right hand side

        self.fig = plt.figure(constrained_layout=True)
        gs = self.fig.add_gridspec(5, 5)
        self.ax_im = self.fig.add_subplot(gs[:, 0:4])
        self.ax_ca = self.fig.add_subplot(gs[0, -1])
        self.ax_force = self.fig.add_subplot(gs[1, -1])
        self.ax_length = self.fig.add_subplot(gs[2, -1])
        self.ax_thin = self.fig.add_subplot(gs[3, -1])
        self.ax_thick = self.fig.add_subplot(gs[4, -1])

    def display_snapshot(self, render_file_string):
        """ Shows the snapshot image file """

        self.ax_im.imshow(plt.imread(render_file_string))
        self.ax_im.axis('off')
        
    def save_fig_as_image(self, output_image_file_string, dpi=300):
        """ Saves current figure as image_file_string """
        
        # Save if required
        if output_image_file_string:
            print('Saving figure to %s' % output_image_file_string)
            # Check path exists
            folder = os.path.dirname(output_image_file_string)
            if not os.path.exists(folder):
                os.makedirs(folder)
            self.fig.savefig(output_image_file_string, dpi=dpi)