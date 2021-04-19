"""@package FiberSim_analysis
Documentation for the FiberSim_analysis module

More stuff here
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from package.modules.half_sarcomere import half_sarcomere as hs

import os
from natsort import natsorted



class analysis():
    """Documentation for the analysis class """

#     More details
      
    def __init__(self):
        """Constructor"""
        
    def animate_cb_distributions(self, hs_status_directory,
                                output_gif_file,
                                frames=[]):
        """Creates an animation showing cb distributions over time """
        import os
        from natsort import natsorted
        import imageio

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
        if (frames):
            hs_files = hs_files[frames[0]:frames[1]]

        # Create an animated gig
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
        
    def animate_filament_lengths(self, hs_status_directory,
                                 output_gif_file,
                                 frames=[]):
        """ Creates an animation showing filament lengths over time """
        import imageio

        print(hs_status_directory)
        print(output_gif_file)

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
        if (frames):
            hs_files = hs_files[frames[0]:frames[1]]
        
        # Create an animated gig
        with imageio.get_writer(output_gif_file, mode='I') as writer:
            # Now loop through files
            for i,f in enumerate(hs_files):
                # Load the json file as an hs object
                h = hs.half_sarcomere(f)
                # Make new image file
                pre, ext = os.path.splitext(f)
                image_file_string = pre + '.png'
                # Draw the filament lengths
                h.draw_filament_lengths(image_file_string)
                del h
                # Add in to the gif
                image = imageio.imread(image_file_string, format='png')
                writer.append_data(image)

    def animate_filament_states(self, hs_status_directory,
                                 output_gif_file,
                                 frames=[]):
        """ Creates an animation showing filament lengths over time """
        import imageio

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
        if (frames):
            hs_files = hs_files[frames[0]:frames[1]]

        hs_files = hs_files[0:1900:5]

        # Create an animated gif
        with imageio.get_writer(output_gif_file, mode='I') as writer:
            # Loop through files
            for i,f in enumerate(hs_files):
                h = hs.half_sarcomere(f)
                # Make new image file
                pre, ext = os.path.splitext(f)
                image_file_string = pre +'.png'
                # Draw the thick filament states
                h.draw_myofilament_states(image_file_string)
                del h
                # Add in to the gif
                image = imageio.imread(image_file_string, format='png')
                writer.append_data(image)


    def display_pCa_data(self, data_directory, y_data_type = "force",
                         output_file_string=""):
        """Pulls x-pCa data from a folder"""
        
        # Get the data files
        data_files = []
        for dir_path,_,filenames in os.walk(data_directory):
            for f in filenames:
                f = os.path.abspath(os.path.join(data_directory, f))
                if (f.endswith(".txt")):
                    data_files.append(f)
        
        # Now sort them
        data_files = natsorted(data_files)

        # Create holders for the data
        pCa = []
        y = []
        
        # Loop through files
        for i,f in enumerate(data_files):
            fs_data = pd.read_table(f)
            pCa.append(fs_data['pCa'].iloc[-1])
            y.append(fs_data[y_data_type].iloc[-1])

        # Fit the data
        from .utilities import fit_pCa_data
        fit_result = fit_pCa_data(pCa, y)

        # Make the figure
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([3,3])
        spec2 = gridspec.GridSpec(nrows = 1, ncols = 1,
                                  figure = f)

        ax1 = f.add_subplot(spec2[0, 0]);
        ax1.plot(pCa, y, 'go')

        # Add in fit
        ax1.plot(fit_result['x_fit'],fit_result['y_fit'],'r-',
                 label='pCa50 = %.3f\nn_H = %.2f' %
                     (fit_result['pCa_50'],fit_result['n_H']))
        ax1.legend()
        ax1.set_xlabel('pCa')
        ax1.set_ylabel(y_data_type)

        if (output_file_string):
            from .utilities import save_figure_to_file
            save_figure_to_file(f, output_file_string)


    def display_data(self, data_file_string,
                     t_limits=[],
                     output_file_string=""):
        """Summary display"""
        fs_data = pd.read_table(data_file_string)
        
        if t_limits:
            t = fs_data['time']
            vi = np.nonzero((t>t_limits[0])&(t<=t_limits[1]))
            fs_data = fs_data.iloc[vi]
        
        no_of_rows = 6
        no_of_cols = 1
        
        f = plt.figure(constrained_layout=True)
        f.set_size_inches([4,8])
        spec2 = gridspec.GridSpec(nrows = no_of_rows, ncols = no_of_cols,
                                  figure = f)
        
        ax1 = f.add_subplot(spec2[0, 0])
        ax1.plot('time', 'force', data=fs_data)
        ax1.set_ylabel('Force (N m^{-2})')
        
        ax2 = f.add_subplot(spec2[1, 0])
        ax2.plot('time', 'hs_length', data=fs_data)
        ax2.set_ylabel('hs length (nm)')
        
        ax3 = f.add_subplot(spec2[2, 0])
        ax3.plot('time', 'a_fil_length', data=fs_data)
        ax3.plot('time', 'm_fil_length', data=fs_data)
        ax3.set_ylabel('Filament\nlengths\n(nm)')
        
        ax4 = f.add_subplot(spec2[3, 0])
        ax4.plot('time', 'a_pop_0', data=fs_data)
        ax4.plot('time', 'a_pop_1', data=fs_data)
        ax4.set_ylabel('a_pops')

        ax5 = f.add_subplot(spec2[4, 0])
        m_pop_labels = [col for col in fs_data.columns if 'm_pop' in col]
        for m in m_pop_labels:
            ax5.plot('time', m, data=fs_data)
        ax5.set_ylabel('m_pops')

        ax6 = f.add_subplot(spec2[5, 0])
        ax6.plot('time', 'c_pop_0', data=fs_data)
        ax6.plot('time', 'c_pop_1', data=fs_data)
        ax6.set_ylabel('c_pops')

        if (output_file_string):
            from .utilities import save_figure_to_file
            save_figure_to_file(f, output_file_string)
        else:
            plt.show()
