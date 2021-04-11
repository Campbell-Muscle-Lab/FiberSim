# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 22:31:21 2021

@author: ken
"""

import os
import re
import json
import time

import pandas as pd

import cv2


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class movie_with_data():
    """ Class to create a movie that includes screen-shots and plots """

    def __init__(self, render_file_string):
        """ Constructor """
        
        render_file_string = 'c:/ken/github/campbellmusclelab/models/fibersim/demo_files/visualization/movie_with_data/movie.json'

        # Hold the file string
        self.render_file_string = render_file_string
        self.parent_path = os.path.dirname(self.render_file_string)


        # Load the data
        with open(render_file_string, 'r') as f:
            self.render_dict = json.load(f)
        print(self.render_dict)

        # Load the results data
        results_file = os.path.join(self.parent_path,
                                    self.render_dict['results']['results_file'])
        
        self.d = pd.read_csv(results_file, sep='\t')
        
        self.deduce_frames()
        
        self.movie_frames =[]
        self.loop_through_frames()
        
        self.write_movie()
    
    def write_movie(self):

        img_array = []
        print('ken was here')
        print(self.movie_frames)
        for f in self.movie_frames:
            print(f)
            im = cv2.imread(f)
            height, width, layers = im.shape
            size = (width,height)
            img_array.append(im)
        
        out = cv2.VideoWriter('c:/temp/ken.mp4', cv2.VideoWriter_fourcc(*'MP4V'), 10, size)
        
        for i in range(len(img_array)):
            out.write(img_array[i])
        out.release()
        print('Written movie')
        

    def deduce_frames(self):
        # Deduces the frames

        # Find the frames
        render_data = self.render_dict['render_batch']['render_jobs'][0]
        print(render_data)
        
        if (render_data['relative_to'] == 'this_file'):
            parent_path = os.path.dirname(self.render_file_string)
            frames_file = os.path.join(parent_path,
                                       render_data['frames_file'])
        
        # Load frames data
        with open(frames_file, 'r') as f:
            frames_dict = json.load(f)
        
        self.top_files = []
        self.top_indices = []
        self.top_c_x = []
        self.bot_files = []
        self.bot_indices = []
        self.bot_c_x = []
        for f in frames_dict['frames']:
            fs = f['image_file']
            vi = re.findall(r'[0-9]+', fs)
            cx = f['camera']['location']['x']
            
            if ('c_zone' in f['image_file']):
                self.top_files.append(f['image_file'])
                self.top_indices.append(int(vi[0]))
                self.top_c_x.append(cx)
            else:
                self.bot_files.append(f['image_file'])
                self.bot_indices.append(int(vi[0]))
                self.bot_c_x.append(cx)

        # Sort the top and bot
        zipped = zip(self.top_indices,
                     self.top_files, self.top_c_x)
        s = sorted(zipped)
        t = zip(*s)
        self.top_indices, self.top_files, self.top_c_x = \
            [list(t1) for t1 in t]

        zipped = zip(self.bot_indices,
                     self.bot_files, self.bot_c_x)
        s = sorted(zipped)
        t = zip(*s)
        self.bot_indices, self.bot_files, self.bot_c_x = \
            [list(t1) for t1 in t]

    def loop_through_frames(self):
        """ Loops through all the frames saving an image for each one """

        for i,fi in enumerate(self.top_indices):
            # if not (i%100 == 0):
            #     continue
            
            # Generate top image
            top_file = os.path.join(self.parent_path,
                                   'renders/hs_c_zone_%i.png' % fi)
            bot_file = os.path.join(self.parent_path,
                                   'renders/hs_i_band_%i.png' % fi)
            
            # self.create_figure()
            
            # self.display_snapshots(top_file, bot_file)
            # self.add_plots(fi)
            # self.add_cameras(fi, i)
            
            temp_file = os.path.join(self.parent_path,
                                     'temp',
                                     'f_%i.png' % fi)
            
            self.movie_frames.append(temp_file)
            
            # self.save_fig_as_image(temp_file)
            
            # plt.close(self.fig)

    def add_plots(self, fi, max_t = 1):
        """ Adds the data """
        
        t = self.d['time'].iloc[0:fi].to_numpy()
        pCa = self.d['pCa'].iloc[0:fi].to_numpy()
        hsl = self.d['hs_length'].iloc[0:fi].to_numpy()
        force = 0.001 * self.d['force'].iloc[0:fi].to_numpy()
        n_off = self.d['a_pop_0'].iloc[0:fi].to_numpy()
        n_on = self.d['a_pop_1'].iloc[0:fi].to_numpy()
        m_srx = self.d['m_pop_0'].iloc[0:fi].to_numpy()
        m_drx = self.d['m_pop_1'].iloc[0:fi].to_numpy()
        m_bound = self.d['m_pop_2'].iloc[0:fi].to_numpy()

        self.ax_pCa.plot(t, pCa)
        self.ax_pCa.set_xlim([0, max_t])
        self.ax_pCa.set_ylim([10, 5])
        self.ax_pCa.set_ylabel('pCa',
                                  rotation=0,
                                  loc='center',
                                  verticalalignment='center',
                                  labelpad = 30)

        self.ax_length.plot(t, hsl)
        self.ax_length.set_xlim([0, max_t])
        self.ax_length.set_ylim([900, 1300])
        self.ax_length.set_ylabel('Half-\nsarcomere\nlength\n(nm)',
                                  rotation=0,
                                  loc='center',
                                  verticalalignment='center',
                                  labelpad = 30)
        
        self.ax_force.plot(t, force)
        self.ax_force.set_xlim([0, max_t])
        self.ax_force.set_ylim([0, 70])
        self.ax_force.set_ylabel('Stress\n(kN m^{-2})',
                                  rotation=0,
                                  loc='center',
                                  verticalalignment='center',
                                  labelpad = 30)

        self.ax_thin.plot(t, n_on, label='N_on')
        self.ax_thin.plot(t, n_off, label='N_off')
        self.ax_thin.set_xlim([0, max_t])
        self.ax_thin.set_ylim([0, 1])
        self.ax_thin.set_ylabel('Thin\nfilament',
                                  rotation=0,
                                  loc='center',
                                  verticalalignment='center',
                                  labelpad = 30)
        self.ax_thin.legend(loc='upper left')

        self.ax_thick.plot(t, m_srx, label='M_SRX')
        self.ax_thick.plot(t, m_drx, label='M_DRX')
        self.ax_thick.plot(t, m_bound, label='M_bound')
        self.ax_thick.set_xlim([0, max_t])
        self.ax_thick.set_ylim([0, 1])
        self.ax_thick.set_ylabel('Thick\nfilament',
                                  rotation=0,
                                  loc='center',
                                  verticalalignment='center',
                                  labelpad = 30)
        self.ax_thick.legend(loc='upper left')
        

    def add_cameras(self, fi, ci):
        # Adds camera
        
        a_end = self.d['a_fil_length'].iloc[fi]
        hsl = self.d['hs_length'].iloc[fi]
        m_length = self.d['m_fil_length'].iloc[fi]
        
        self.ax_cam_top.plot([0, a_end],[0.25, 0.25], 'k-',
                             linewidth=3)
        self.ax_cam_top.plot([hsl-m_length, hsl],[1, 1], 'k-',
                             linewidth=10)
        self.ax_cam_top.plot([0,0],[0,2],'k-',
                             linewidth=3)
        self.ax_cam_top.plot([hsl,hsl],[0,2],'k-',
                             linewidth=10)
        self.ax_cam_top.set_xlim([0, 1300])
        self.ax_cam_top.set_ylim([0, 2])
        self.ax_cam_top.text(self.top_c_x[ci], 1.5, 'Cam->')
        self.ax_cam_top.axis('off')

        self.ax_cam_bot.plot([0, a_end],[0.25, 0.25], 'k-',
                             linewidth=3)
        self.ax_cam_bot.plot([hsl-m_length, hsl],[1, 1], 'k-',
                             linewidth=10)
        self.ax_cam_bot.plot([0,0],[0,2],'k-',
                             linewidth=3)
        self.ax_cam_bot.plot([hsl,hsl],[0,2],'k-',
                             linewidth=10)
        self.ax_cam_bot.set_xlim([0, 1300])
        self.ax_cam_bot.set_ylim([0, 2])
        self.ax_cam_bot.text(self.bot_c_x[ci], 1.5, 'Cam->')
        
        self.ax_cam_bot.axis('off')


    def create_figure(self):
        """ Initialises figure """

        # Make a figure iwth 1 large panel to show images and
        # a column of data panels on right hand side
        
        self.fig = plt.figure(constrained_layout=False)
        self.fig.set_size_inches([9, 9])
        self.gs = self.fig.add_gridspec(10, 10)
        self.ax_im_1 = self.fig.add_subplot(self.gs[1:5, 0:6])
        self.ax_im_2 = self.fig.add_subplot(self.gs[6:11, 0:6])
        
        self.ax_cam_top = self.fig.add_subplot(self.gs[0, 0:6])
        self.ax_cam_bot = self.fig.add_subplot(self.gs[5, 0:6])
        
        self.ax_pCa = self.fig.add_subplot(self.gs[0:2, 7:11])
        self.ax_length = self.fig.add_subplot(self.gs[2:4, 7:11])
        self.ax_force = self.fig.add_subplot(self.gs[4:6, 7:111])
        self.ax_thin = self.fig.add_subplot(self.gs[6:8, 7:11])
        self.ax_thick = self.fig.add_subplot(self.gs[8:10, 7:11])

        self.gs.tight_layout(self.fig)


    def display_snapshots(self, top_image_file, bot_image_file):
        """ Shows the snapshot image file """

        self.ax_im_1.imshow(plt.imread(top_image_file))
        self.ax_im_1.axis('off')
        
        self.ax_im_2.imshow(plt.imread(bot_image_file))
        self.ax_im_2.axis('off')
        
    def save_fig_as_image(self, output_image_file_string, dpi=100):
        """ Saves current figure as image_file_string """
        
        # Save if required
        if output_image_file_string:
            print('Saving figure to %s' % output_image_file_string)
            # Check path exists
            folder = os.path.dirname(output_image_file_string)
            if not os.path.exists(folder):
                os.makedirs(folder)
            self.fig.savefig(output_image_file_string, dpi=dpi)
            
if __name__ == "__main__":
    c = movie_with_data('a')