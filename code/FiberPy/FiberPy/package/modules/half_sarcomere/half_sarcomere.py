# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:11:28 2020
Modified on Monday Jun 22 by Sarah

@author: kscamp3
"""

import os

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.append("..") 

from thin_filament import thin_filament as thin
from thick_filament import thick_filament as thick


class half_sarcomere():
    """Half-sarcomere class"""
    
    def __init__(self, json_file_string):
       
        with open(json_file_string) as json_file:
           json_data = json.load(json_file)
 
           # Load hs data
           hs_data = json_data['hs_data']
           self.hs_id = hs_data['hs_id']
           self.time = hs_data['time']
           self.hs_length = hs_data['hs_length']
           self.hs_force = hs_data['hs_force']
           self.pCa = hs_data['pCa']
           self.m_nodes_per_thick_filament = \
               hs_data['m_nodes_per_thick_filament']
           self.a_nodes_per_thin_filament = \
               hs_data['a_nodes_per_thin_filament']
               
           self.t_attach_a_node = hs_data['titin']['t_attach_a_node']
           self.t_attach_m_node = hs_data['titin']['t_attach_m_node']
           self.t_slack_length = hs_data['titin']['t_slack_length']
           self.t_k_stiff = hs_data['titin']['t_k_stiff']
           self.cb_extensions = hs_data['cb_extensions']

           #Load thick filaments
           thick_fil_data = hs_data['thick']
           self.thick_fil = []
           for t in thick_fil_data:
               self.thick_fil.append(thick.thick_filament(t))

           #Load thin filaments
           thin_fil_data = hs_data['thin']
           self.thin_fil = []
           for t in thin_fil_data:
               self.thin_fil.append(thin.thin_filament(t))

    def draw_cb_distributions(self, output_file_string=""):
        """ Draws cb_x histograms """

        # Loops through thick filaments looking for bound heads
        cb_displacements = []
        for t in self.thick_fil:
            for i, bound_f in enumerate(t.cb_bound_to_a_f):
                if (bound_f >= 0):
                    nf = t.cb_bound_to_a_f[i]
                    nn = t.cb_bound_to_a_n[i]
                    print("nf: %i nn %i" % (nf, nn))
                    x1 = t.cb_x[i]
                    x2 = self.thin_fil[bound_f].bs_x[t.cb_bound_to_a_n[i]]
                    cb_displacements.append(x1-x2);

        # Generate the histogram
        b = np.linspace(-15.5, 15.5, num=32)
        y, bin_edges = np.histogram(cb_displacements, bins=b)

        # Set up for display
        x = b[0:-1]+0.5
        y = y / (len(self.thick_fil)*self.thick_fil[0].m_no_of_cbs)

        # Display
        no_of_rows = 2
        no_of_cols = 2

        f = plt.figure(constrained_layout=True)
        f.set_size_inches([4,4])
        spec = gridspec.GridSpec(nrows=no_of_rows, ncols=no_of_cols,
                                 figure=f)

        ax1 = f.add_subplot(spec[0,0])
        ax1.plot(x,y);
        ax1.set_ylim(0,0.15)
        ax1.set_xlabel("Cross-bridge lenth (nm)")
        ax1.set_ylabel("Proportion\nof\cross-bridges")

        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()

    def draw_filaments(self, output_file_string=""):
        
        no_of_rows = 2;
        no_of_cols = 2;
        
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([4,4])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)

        y_m = np.ones((self.thick_fil[0].m_no_of_cbs, 1))
        y_c = np.ones((self.thick_fil[0].c_no_of_pcs, 1))
        y_a = 1+np.ones((self.thin_fil[0].a_no_of_bs, 1))
        
        ax1 = f.add_subplot(spec[0,0]);
        ax1.plot(self.thick_fil[0].cb_x, y_m, 'b+')
        ax1.plot(self.thick_fil[0].pc_x, y_c, 'ro')
        ax1.plot(self.thin_fil[0].bs_x, y_a, 'g+')

    def draw_myofilament_states(self, output_file_string=""):

        # Map thick filaments
        for i, t in enumerate(self.thick_fil):
            if (i==0):
                m = np.zeros([3,len(self.thick_fil),len(t.cb_state)])
                m_x = t.cb_x
            for s in range(0,3):
                m[s,i,:] = (t.cb_state/(s+1)) * (t.cb_state == (s+1))

        # Map thin filaments
        for i, t in enumerate(self.thin_fil):
            if (i==0):
                a = np.zeros([2,len(self.thin_fil),len(t.bs_state)])
                a_x = t.bs_x
            for s in range(0,2):
                if (s==0):
                    vi = np.nonzero(t.bs_state==0)
                    a[s,i,vi]=1
                else:
                    a[s,i,:] = (t.bs_state/(s)) * (t.bs_state == (s))

        no_of_rows = 1
        no_of_cols = 1
        
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([4,2])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)
        
        c = [[1, 0, 0],[0, 1, 0],[0, 0, 1],[0.5, 0.5, 0],[0.5, 0, 0.5]]
        labels = ['M_OFF','M_ON','M_attached','A_OFF','A_ON']

        ax_0_0 = f.add_subplot(spec[0, 0])
        for s in range(0,3):
            ax_0_0.plot(m_x, 2+np.mean(m[s,:,:],0),color=c[s],label=labels[s])
        for s in range(0,2):
            ax_0_0.plot(a_x, np.mean(a[s,:,:],0), color=c[s+3],label=labels[s+3])
        ax_0_0.legend()
        ax_0_0.set_ylim(0,3)
        
        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()


    def draw_thick_filament_states(self, output_file_string=""):
        
        no_of_rows = 3
        no_of_cols = 1

        f = plt.figure(constrained_layout=True)
        f.set_size_inches([8,4])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)

        for i, t in enumerate(self.thick_fil):
            if (i==0):
                s1 = np.zeros(len(t.cb_state))
                s2 = s1
                s3 = s1
            s1 = s1 + (t.cb_state * (t.cb_state == 1))
            s2 = s2 + (t.cb_state * (t.cb_state == 2))
            s3 = s3 + (t.cb_state * (t.cb_state == 3))
        
        s1 = s1 / (len(self.thick_fil))
        s2 = s2 / (2 * len(self.thick_fil))
        s3 = s3 / (3 * len(self.thick_fil))
        
        ax_0_0 = f.add_subplot(spec[0, 0])
        ax_0_0.plot(s1,'b-')
        ax_0_0.set_ylim(0,1)
        ax_0_0.set_ylabel('OFF state')

        ax_1_0 = f.add_subplot(spec[1, 0])
        ax_1_0.plot(s2,'b-')
        ax_1_0.set_ylim(0,1)
        ax_1_0.set_ylabel('ON state')

        ax_2_0 = f.add_subplot(spec[2, 0])
        ax_2_0.plot(s3,'b-')
        ax_2_0.set_ylim(0,1)
        ax_2_0.set_ylabel('Attached state')
        ax_2_0.set_xlabel('Cross-bridge position on thick filament')
        
        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()


    def draw_filament_lengths(self, output_file_string=""):
        """ Draws cb_x and bs_x positions """
        
        no_of_rows = 3;
        no_of_cols = 2;
        
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([8,4])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)

        # Thick filaments
        ax2 = f.add_subplot(spec[1,0])
        for t in self.thick_fil:
            ax2.plot(t.cb_x)
        ax2.set_ylabel('Cross-bridge\npositions\n(nm)',
                       rotation=90)
        ax2.set_ylim(0,1800)

        ax4 = f.add_subplot(spec[2,0])
        for t in self.thick_fil:
            ax4.plot(-np.diff(t.cb_x[0:-1:6]))
        ax4.set_ylabel('Inter-\ncrown\ndistance\n(nm)')
        # ax4.set_ylim(14, 18)

        # Thin filaments
        ax3 = f.add_subplot(spec[1,1])
        for t in self.thin_fil:
            ax3.plot(t.bs_x)
        ax3.set_ylabel('Binding site\npositions\n(nm)',
                       rotation=90)

        ax5 = f.add_subplot(spec[2,1])
        for t in self.thin_fil:
            ax5.plot(np.diff(t.bs_x[0:-1:2]))
        ax5.set_ylabel('Inter-site\ndistance\n(nm)',
                       rotation=90)
        # ax5.set_ylim(5, 6)

        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()
        
        
    def save_figure_to_file(self, f, im_file_string, dpi=150, verbose=1):
        """ Writes an image to file """
        
        dir_path = os.path.dirname(im_file_string)
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
            
        if (verbose):
            print('Saving figure to: %s' % im_file_string)
            
        f.savefig(im_file_string, dpi=dpi)