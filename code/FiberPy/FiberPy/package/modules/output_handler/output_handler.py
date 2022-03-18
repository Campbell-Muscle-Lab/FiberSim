# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 17:17:59 2020

@author: ken
"""

import json
import os

import numpy as np
import pandas as pd

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from ..display.multi_panel import multi_panel_from_flat_data


class output_handler():
    """ Class for handling simulation output """

    def __init__(self, output_handler_file_string,
                 sim_data=[],
                 sim_results_file_string=[]):

        # Display
        print("Initialising FiberPy output_handler")

        # Check for output_handler file
        if (output_handler_file_string == []):
            print('No output handler file specified. Cannot write output')
            return

        # Load the output handler structure as a dict
        with open(output_handler_file_string, 'r') as f:
            self.oh_data = json.load(f)

        # Check for a sim_results file, overwriting sim data that may
        # have been passed in
        if sim_results_file_string:
            print('Loading sim data from %s' % sim_results_file_string)
            sim_data = pd.read_csv(sim_results_file_string,
                                   delimiter='\t')

        # Check we have data to do something with
        if not isinstance(sim_data, pd.DataFrame):
            print('No simulation data available')
            return

        # User defined files
        if ('templated_images' in self.oh_data):
            if ('relative_to' in self.oh_data['templated_images'][0]):
                user_defined = self.oh_data['templated_images']
                for ud in user_defined:
                    if (not ud['relative_to']):
                        template_fs = os.path.abspath(ud['template_file_string'])
                        output_fs = os.path.abspath(ud['output_file_string'])
                    elif (ud['relative_to'] == 'this_file'):
                        base_directory = \
                            Path(output_handler_file_string).parent.absolute()
                        template_fs = os.path.join(base_directory,
                                                   ud['template_file_string'])
                        output_fs = os.path.join(base_directory,
                                                 ud['output_file_string'])
                    else:
                        base_directory = ud['relative_to']
                        template_fs = os.path.join(base_directory,
                                                   ud['template_file_string'])
                        output_fs = os.path.join(base_directory,
                                                 ud['output_file_string'])
    
            if (not 'output_image_formats' in ud):
                ud['output_image_formats'] = ['png']

            self.create_image_from_template(
            sim_data,
            template_fs,
            output_fs,
            ud['output_image_formats'])

    def create_image_from_template(self,
                                   sim_data,
                                   template_file_string,
                                   output_image_file_string,
                                   image_formats):

        if not os.path.isabs(template_file_string):
            template_file_string = os.path.join(os.getcwd(),
                                                template_file_string)
        if not os.path.isabs(output_image_file_string):
            output_image_file_string = os.path.join(os.getcwd(),
                                                    output_image_file_string)
        fig, ax = multi_panel_from_flat_data(
                        pandas_data=sim_data,
                        template_file_string=template_file_string,
                        output_image_file_string=[])
        
        for f in image_formats:
            ofs = ('%s.%s' % (output_image_file_string, f))
            print('Writing image to: %s' % ofs)
            fig.savefig(ofs, dpi=200, bbox_inches='tight')

        plt.close(fig)

    def animate_cb_distributions(self,
                                 cb_dump_file_string=[],
                                 output_image_file_string=[],
                                 skip_frames=1):
        """ Animates a cb distribution """

        import imageio

        # Open file and pull x values
        with open(cb_dump_file_string, 'r') as f:
            x_strings = f.readline().split('\t')
        f.close()
        x = []
        for i, xs in enumerate(x_strings[1:]):
            x.append(float(xs[1:]))
        # Now load distribs as numpy arrray
        cb_dump = np.loadtxt(cb_dump_file_string, skiprows=1)
        cb_distribs = cb_dump[:, 1:]
        t = cb_dump[:, 0]
        max_pop = np.amax(cb_distribs)

        temp_image_file_string = 'temp.png'
        print('Animating cross-bridge distribution')
        with imageio.get_writer(output_image_file_string, mode='I') \
                as writer:
            for i in np.arange(0, np.shape(cb_distribs)[0],
                               skip_frames):
                print(('Frame: %.0f' % i), end=' ', flush=True)
                self.draw_cb_distribution(x, cb_distribs[i, :],
                                          t[i], 1.2*max_pop,
                                          temp_image_file_string)
                image = imageio.imread(temp_image_file_string, format='png')
                writer.append_data(image)
            print('Animation built')
            print('Animation written to %s' % output_image_file_string)
        os.remove(temp_image_file_string)
