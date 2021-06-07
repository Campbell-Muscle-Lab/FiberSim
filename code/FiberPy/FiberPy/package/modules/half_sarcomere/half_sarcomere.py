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

class half_sarcomere(dict):
    """Half-sarcomere class"""
    
    def __init__(self, json_file_string):

        with open(json_file_string) as json_file:
           json_data = json.load(json_file)
 
           self['hs_data'] = json_data['hs_data']
           self['titin'] = json_data['titin']
           self['thick'] = json_data['thick']
           self['thin'] = json_data['thin']

    def draw_cb_distributions(self, output_file_string=""):
        """ Draws cb_x histograms 
        
        Parameters
        ----------
        output_file_string : string, optional
            Name of the (optional) output figure file. 

        Returns
        -------
        None.
        """

        # Loops through thick filaments looking for bound heads
        cb_displacements = []
        for t in self["thick"]:
            for i, bound_f in enumerate(t['cb_bound_to_a_f']):
                if (bound_f >= 0):
                    nf = t['cb_bound_to_a_f'][i]
                    nn = t['cb_bound_to_a_n'][i]
                    x1 = t['cb_x'][i]
                    x2 = self['thin'][nf]['bs_x'][nn]
                    cb_displacements.append(x1-x2);

        # Generate the histogram
        b = np.linspace(-15, 15, num=200)
        y, bin_edges = np.histogram(cb_displacements, bins=b)

        # Set up for display
        x = b[0:-1]+0.5
        y = y / (len(self["thick"])*self["thick"][0]["m_no_of_cbs"])

        # Display
        no_of_rows = 1
        no_of_cols = 1

        f = plt.figure(constrained_layout=True)
        f.set_size_inches([4,4])
        spec = gridspec.GridSpec(nrows=no_of_rows, ncols=no_of_cols,
                                 figure=f)

        ax1 = f.add_subplot(spec[0,0])
        ax1.plot(x,y);
        #ax1.set_ylim(0,0.15)
        ax1.set_xlabel("Cross-bridge stretch (nm)")
        ax1.set_ylim([0, 0.2])
        ax1.set_ylabel("Proportion\nof\ncross-bridges")

        plt.show()
        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()
    
    def draw_filaments(self, output_file_string=""):
        """
        Draw the first thick and thin filaments

        Parameters
        ----------
        output_file_string : string, optional
            Name of the (optional) output figure file. 

        Returns
        -------
        None.

        """
        
        no_of_rows = 1;
        no_of_cols = 1;
        
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([4,4])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)

        y_m = np.ones((self["thick"][0]["m_no_of_cbs"], 1))
        y_c = np.ones((self["thick"][0]["c_no_of_pcs"], 1))
        y_a = 1+np.ones((self["thin"][0]["a_no_of_bs"], 1))
        
        ax1 = f.add_subplot(spec[0,0]);
        ax1.plot(self["thick"][0]["cb_x"], y_m, 'b+', label = "thick")
        #ax1.plot(self["thick"][0]["pc_x"], y_c, 'ro', label = "C-protein")
        ax1.plot(self["thin"][0]["bs_x"], y_a, 'g+', label = "thin")
        ax1.set_yticks([])
        ax1.set_xlabel("x axis (nm)")
        ax1.set_ylim([0,2.1])
        
        ax1.legend(loc = "best")
        
               
    def draw_myofilament_states(self, no_of_states = 2, output_file_string=""):
        """
        Plots the actin and myosin states as a function of their x-positions
        
        Parameters
        ----------
        no_of_states : int, optional
            Number of myosin states. The default is 2.
        output_file_string : string, optional
            Name of the (optional) output figure file. 

        Returns
        -------
        None.

        """        
    

        # Map thick filaments
        for i, t in enumerate(self["thick"]):
            if (i==0):
                m = np.zeros([no_of_states,len(self["thick"]),len(t["cb_state"])])
                m_x = t["cb_x"]
            for s in range(0,no_of_states):
                m[s,i,:] = [state/(s+1) * (state == (s+1)) for state in t["cb_state"]]

        # Map thin filaments
        for i, t in enumerate(self["thin"]):
            if (i==0):
                a = np.zeros([2,len(self["thin"]),len(t["bs_state"])])
                a_x = t["bs_x"]
            for s in range(0,2):
                if (s==0):
                    vi = np.nonzero([x==1 for x in t["bs_state"]])
                    a[s,i,vi]=1
                else:
                    a[s,i,:] = [ (x/(s + 1)) * (x == (s + 1)) for x in t["bs_state"] ]

        no_of_rows = 1
        no_of_cols = 1
        
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([4,2])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)
        
        c = [[1, 0, 0],[0, 1, 0],[0, 0, 1],[1, 0.5, 0.5],[0.5, 0, 0.5]]
        
        labels = []
        for s in range(0, no_of_states):
            labels.append(f"M_{s+1}")
        labels.append('A_OFF')
        labels.append('A_ON')

        ax_0_0 = f.add_subplot(spec[0, 0])
        ax_0_0.set_title("Myofilament states")
        for s in range(0,no_of_states):
            ax_0_0.plot(m_x, no_of_states + np.mean(m[s,:,:],0),color=c[s],label=labels[s])
        for s in range(0,2):
            ax_0_0.plot(a_x, np.mean(a[s,:,:],0), color=c[s+no_of_states],label=labels[s+no_of_states])
        ax_0_0.legend()
        ax_0_0.set_ylim(0,no_of_states+1)
        ax_0_0.set_xlabel("x position (nm)")
        
        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()


    def draw_thick_filament_states(self, no_of_states = 2, output_file_string=""):
        """
        Plots the myosin states 

        Parameters
        ----------
        no_of_states : int, optional
            Number of myosin states. The default is 2.
        output_file_string : string, optional
            Name of the (optional) output figure file. 

        Returns
        -------
        None.

        """
        
        no_of_rows = no_of_states
        no_of_cols = 1

        f = plt.figure(constrained_layout=True)
        f.set_size_inches([8,4])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)    
                
        for i, t in enumerate(self["thick"]):
            
            if (i==0):
                s = np.zeros([len(t["cb_state"]), no_of_states])    
                
            for state in range (0, no_of_states):
                s[:, state] = s[:, state] + [x * (x == state + 1) for x in t["cb_state"]]
        
        for state in range (0, no_of_states):
            
            s[:, state] = s[:, state]/( (state + 1) * len(self["thick"])) 
                
            ax = f.add_subplot(spec[state, 0])
            ax.plot(s[:, state],'-')
            ax.set_ylim(0,1)
            ax.set_ylabel(f'State {state+1}')
            
            if state == no_of_states - 1:
                ax.set_xlabel('# Cross-bridge')
        f.suptitle("Myosin states")
                
        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()
            
    def draw_c_protein_states(self, no_of_states = 2, output_file_string=""):
        """
        Plots the C-protein states

        Parameters
        ----------
        no_of_states : int, optional
            Number of C-protein states. The default is 2.
        output_file_string : string, optional
            Name of the (optional) output figure file. 

        Returns
        -------
        None.

        """
                
        no_of_rows = no_of_states
        no_of_cols = 1

        f = plt.figure(constrained_layout=True)
        f.set_size_inches([8,4])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)    
                
        for i, t in enumerate(self["thick"]):
            
            if (i==0):
                s = np.zeros([len(t["pc_state"]), no_of_states])    
                
            for state in range (0, no_of_states):
                s[:, state] = s[:, state] + [x * (x == state + 1) for x in t["pc_state"]]
        
        for state in range (0, no_of_states):
            
            s[:, state] = s[:, state]/( (state + 1) * len(self["thick"])) 
                
            ax = f.add_subplot(spec[state, 0])
            ax.plot(s[:, state],'-')
            ax.set_ylim(0,1.1)
            ax.set_ylabel(f'State {state +1}')
            
            if state == no_of_states - 1:
                ax.set_xlabel('# C-protein')
                
        f.suptitle("C-protein states")
                
        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()
            
    def draw_filament_lengths(self, output_file_string=""):
        
        """
        Draws cb_x and bs_x positions

        Parameters
        ----------
        output_file_string : string, optional
            Name of the (optional) output figure file. 

        Returns
        -------
        None.

        """
        
        no_of_rows = 2;
        no_of_cols = 2;
        
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([8,5])
        spec = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                 figure = f)
        
        """
        # Node force 
        
        ax1 = f.add_subplot(spec[0,0])
        for t in self.thick_fil:
            ax1.plot(t.node_forces)         
        ax1.set_yticks([0.0, 0.1])
        ax1.set_yticklabels(['1','2'])
        ax1.set_ylabel('Node force',
                       rotation=90)
        """

        # Thick filaments
        ax2 = f.add_subplot(spec[0,0])
        for t in self["thick"]:
            ax2.plot(t["cb_x"])
        ax2.set_ylabel('Cross-bridge\npositions\n(nm)',
                       rotation=90)
        ax2.set_ylim(0,1800)
        ax2.set_xlabel("# cb")
        
        thin_0 = self["thin"][0]
        thick_0 = self["thick"][0] 
        
        # Beginning of overlap for the thick filament
                      
        overlap_start = max(thin_0["bs_x"])    
        cb_overlap_start = min(thick_0["cb_x"], key=lambda x:abs(x-overlap_start))  
        thick_node = np.where([x == cb_overlap_start for x in thick_0["cb_x"]])[0]/6
        
        ax4 = f.add_subplot(spec[1,0])
        for t in self["thick"]:
            ax4.plot(-np.diff(t["cb_x"][0:-1:6]), label = "_")
            
        ax4.axvline(x=54-thick_node[0], label = "Overlap starts")
        ax4.set_ylabel('Inter-\ncrown\ndistance\n(nm)')
        ax4.set_ylim(13.5,13.51)
        ax4.legend(loc = "upper center")
        ax4.set_xlabel("# thick node")

        # Thin filaments
        ax3 = f.add_subplot(spec[0,1])
        for t in self["thin"]:
            ax3.plot(t["bs_x"])
        ax3.set_ylabel('Binding site\npositions\n(nm)',
                       rotation=90)
        
        ax3.set_xlabel("# bs")
        
        # Beginning of overlap for the thin filament
             
        overlap_start = min(thick_0["cb_x"])        
        bs_overlap_start = min(thin_0["bs_x"], key=lambda x:abs(x-overlap_start))       
        thin_node = np.where([x == bs_overlap_start for x in thin_0["bs_x"]])[0]/2

        ax5 = f.add_subplot(spec[1,1])
        for t in self["thin"]:
            ax5.plot(np.diff(t["bs_x"][0:-1:2]))
            
        ax5.axvline(x=thin_node[0])
        ax5.set_ylabel('Inter-site\ndistance\n(nm)',
                       rotation=90)
        # ax5.set_ylim(5, 6)
        ax5.set_xlabel("# thin node")
        
        f.suptitle("Spatial distributions")

        # Image output
        if (output_file_string):
            self.save_figure_to_file(f, output_file_string)
            plt.close()
      
                    
    def save_figure_to_file(self, f, im_file_string, dpi=150, verbose=1):
        """
        Writes an image to file 

        Parameters
        ----------
        f : figure
            Figure to be saved to a file.
        im_file_string : string
            Name of the figure file.
        dpi : int, optional
            Resolution of the image file. The default is 150.
        verbose : int, optional
            Prints the file path if verbose = 1. The default is 1.

        Returns
        -------
        None.

        """
        
        dir_path = os.path.dirname(im_file_string)
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
            
        if (verbose):
            print('Saving figure to: %s' % im_file_string)
            
        f.savefig(im_file_string, dpi=dpi)

   