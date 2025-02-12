# -*- coding: utf-8 -*-
"""
Created on Wed May 19 16:30:00 2021

@author: kscamp3
"""

import os
import sys
import json
import copy

import numpy as np
import pandas as pd

from pathlib import Path

from collections import defaultdict

from io import StringIO

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from natsort import natsorted

try:
    from package.modules.analysis import curve_fitting as cv
    from package.modules.utilities import utilities as ut
except:
    this_dir = str(Path(os.path.dirname(__file__)).resolve())
   
    sys.path.append(os.path.join(this_dir, '../analysis'))
    import curve_fitting as cv

    sys.path.append(os.path.join(this_dir, '../../'))
    import modules.utilities.utilities as ut

def default_formatting():
    formatting = dict()
    formatting['data_linewidth'] = 1
    formatting['fontname'] = 'Arial'
    formatting['marker_size'] = 8
    formatting['marker_symbols'] = ['o','s','^','v','<','>']
    formatting['fill_styles'] = ['none', 'none', 'none', 'none', 'none', 'none']
    formatting['line_styles'] = ['-','-','-','-','-','-']
    formatting['marker_edge_width'] = 1
    formatting['high_pCa_tick'] = 8.0
    formatting['high_pCa_span'] = 0.2
    formatting['low_pCa_ticks'] = [6.5, 5.5, 4.5]
    formatting['low_pCa_pad'] = 0.2
    formatting['x_label_fontsize'] = 12
    formatting['x_label_pad'] = 15
    formatting['tick_labels_fontsize'] = 12
    formatting['x_scaling_factor'] = 1
    formatting['y_scaling_factor'] = 1
    formatting['x_normalized_to_max'] = False
    formatting['y_normalized_to_max'] = False
    formatting['x_field'] = 'pCa'
    formatting['x_axis_label'] = 'x_axis_label'
    formatting['y_axis_label'] = 'y_axis_label'
    formatting['y_label_pad'] = 15
    formatting['y_label_fontsize'] = 12
    formatting['y_label_rotation'] = 0
    formatting['table_y_anchor'] = 0.95
    formatting['table_y_spacing'] = 0.1
    formatting['table_x_spacing'] = 0.5
    formatting['table_fontsize'] = 11
    formatting['color_set'] = ['tab:blue', 'tab:orange', 'tab:green',
                               'tab:red', 'tab:purple', 'tab:brown',
                               'tab:pink', 'tab:gray', 'tab:olive',
                               'tab:cyan']
    formatting['labels'] = []
    formatting['legend_location'] = 'upper left'
    formatting['legend_bbox_to_anchor'] = [1.05, 1]
    formatting['legend_fontsize'] = 9
    formatting['legend_handlelength'] = 1


    return formatting

def default_layout():
    layout = dict()
    layout['fig_width'] = 3.5
    layout['panel_height'] = 1
    layout['top_margin'] = 0.1
    layout['bottom_margin'] = 0.1
    layout['left_margin'] = 0.1
    layout['right_margin'] = 0.1
    layout['grid_wspace'] = 0.75
    layout['grid_hspace'] = 0.5

    return layout

def plot_y_pCa_data(plot_data, ax=[], y_ticks=[],
                    formatting_overrides=[]):
    """ Plots an array of pCa-y-data """
    
    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if not (formatting_overrides == []):
        for entry in formatting_overrides:
            formatting[entry] = formatting_overrides[entry]
    
    # Make an axes if required
    if (ax == []):
        # Set up current figure
        fig = plt.figure(constrained_layout=False)
        spec = fig.add_gridspec(nrows=1, ncols=1,
                          left=0.3, right=0.95, wspace=0.1,
                          bottom = 0.2)
        fig.set_size_inches([3.5, 2.5])
        
        ax.append(fig.add_suplot(spec[0, 0]))
    
    # Work out how many curves there are
    if (isinstance(plot_data, pd.DataFrame)):
        # There is only one, but we need the data in an array
        no_of_curves = 1
        
        # Make the plot_data, an array of length 1
        pd_temp = []
        pd_temp.append(plot_data)
        
        del plot_data
        plot_data = pd_temp.copy()
    else:
        no_of_curves = len(plot_data)
    
    # Create an empty dictionary
    pd_fit = {}
    
    # Hold min and max values
    min_y = np.inf
    max_y = -np.inf
    
    for i in range(no_of_curves):
        
        # Pull off the values
        pCa = plot_data[i]['pCa']
        y = plot_data[i]['y']
        
        # Calculate the fit
        pd_fit[i] = cv.fit_pCa_data(pCa, y)
        
        x_fit = pd_fit[i]['x_fit']
        y_fit = pd_fit[i]['y_fit']
        
        # Now do some subsitutions
        vi = np.nonzero(pCa >= formatting['high_pCa_tick'])
        pCa[vi] = formatting['high_pCa_tick']
        
        vi = np.nonzero(x_fit >= formatting['high_pCa_tick'])
        x_fit[vi] = formatting['high_pCa_tick']
        
        ax.plot(pCa, y,
                formatting['marker_symbols'][i],
                color = formatting['color_set'][i],
                fillstyle = formatting['fill_styles'][i])
        
        ax.plot(x_fit, y_fit,
                color = formatting['color_set'][i])
            
        max_y = max(max_y, np.amax(y))
        min_y = min(min_y, np.amin(y))
        
    # Set y_ticks
    if (y_ticks == []):
        min_y = ut.multiple_less_than(min_y,
                                      0.05*np.power(10, np.ceil(np.log10(min_y))))
        max_y = ut.multiple_greater_than(max_y,
                                         0.05*np.power(10, np.ceil(np.log10(max_y))))
        y_ticks = [min(min_y, 0), max_y]
    
    # Add in data
    y_anchor = y_ticks[0] + \
        (formatting['table_y_anchor'] * (y_ticks[-1] - y_ticks[0]))
    y_spacing = formatting['table_y_spacing'] * (y_ticks[-1] - y_ticks[0])
    x_anchor = formatting['high_pCa_tick'] - 0.7
    
    ax.text(x_anchor,
            y_anchor,
            'pCa$\\mathregular{_{50}}$',
            fontfamily=formatting['fontname'],
            fontsize=formatting['table_fontsize'],
            horizontalalignment='center',
            verticalalignment='center',
            clip_on=False)
    
    ax.text(x_anchor - (1.5 * formatting['table_x_spacing']),
            y_anchor,
           'n$\\mathregular{_{H}}$',
           fontfamily=formatting['fontname'],
           fontsize=formatting['table_fontsize'],
           horizontalalignment='center',
           verticalalignment='center',
           clip_on=False)
    
    for i in range(no_of_curves):
        y_anchor = y_anchor - y_spacing
        
        ax.plot(x_anchor + 1.25 * formatting['table_x_spacing'],
                y_anchor,
                formatting['marker_symbols'][i],
                color = formatting['color_set'][i],
                fillstyle = formatting['fill_styles'][i])
        
        ax.text(x_anchor,
                y_anchor,
                '%.2f' % pd_fit[i]['pCa_50'],
                fontfamily=formatting['fontname'],
                fontsize=formatting['table_fontsize'],
                horizontalalignment='center',
                verticalalignment='center',
                clip_on=False)
        
        ax.text(x_anchor - (1.5 * formatting['table_x_spacing']),
                y_anchor,
               '%.2f' % pd_fit[i]['n_H'],
               fontfamily=formatting['fontname'],
               fontsize=formatting['table_fontsize'],
               horizontalalignment='center',
               verticalalignment='center',
               clip_on=False)
    
    # Formatting
    ax.invert_xaxis()
    if not (y_ticks == []):
        ax.set_ylim(y_ticks)
    
        
        
    # # Tidy up axis
    # ax_left.set_xlim(formatting['high_pCa_tick'] +
    #                  formatting['high_pCa_span']*np.asarray([0.5, -0.5]))
    # right_lims = [formatting['low_pCa_ticks'][0] + formatting['low_pCa_pad'],
    #               formatting['low_pCa_ticks'][-1]]
    # ax_right.set_xlim(right_lims)
    
    # ax_left.spines['right'].set_visible(False)
    # ax_left.spines['top'].set_visible(False)
    # ax_right.spines['left'].set_visible(False)
    # ax_right.spines['right'].set_visible(False)
    # ax_right.spines['top'].set_visible(False)
    # ax_right.tick_params('y', left=False, labelleft=False)
    
    # # Ticks
    # ax_left.set_xticks([formatting['high_pCa_tick']])
    # ax_left.set_ylabel(formatting['y_axis_label'],
    #                    loc='center',
    #                    verticalalignment='center',
    #                    labelpad=formatting['y_label_pad'],
    #                    fontfamily=formatting['fontname'],
    #                    fontsize=formatting['y_label_fontsize'],
    #                    rotation=formatting['y_label_rotation'])
    
    # ax_right.set_xlabel('pCa',
    #                    loc='center',
    #                    verticalalignment='center',
    #                    labelpad=formatting['x_label_pad'],
    #                    fontfamily=formatting['fontname'],
    #                    fontsize=formatting['x_label_fontsize'])
    # ax_right.set_xticks(formatting['low_pCa_ticks'])
    
    # # Set the y_ticks
    # y_ticks = [0, ut.multiple_greater_than(max_y,
    #                                         0.05*np.power(10, np.ceil(np.log10(max_y))))]
    
    
    # for a in [ax_left, ax_right]:
    #     a.set_xticklabels(a.get_xticks(),
    #                       fontsize=formatting['tick_labels_fontsize'],
    #                       fontfamily=formatting['fontname'])
    #     a.set_ylim(y_ticks)
    #     a.set_yticks(y_ticks)
    #     a.set_yticklabels(a.get_yticks(),
    #                       fontsize=formatting['tick_labels_fontsize'],
    #                       fontfamily=formatting['fontname'])
    
    #     if formatting['y_normalized_to_max']:
    #         y_ticks = [0,1]
    #         a.set_ylim([0,1.03])
    #         a.spines['left'].set_bounds(y_ticks)
    #         a.set_yticks(y_ticks)
    #         a.set_yticklabels(a.get_yticks(),
    #             fontsize=formatting['tick_labels_fontsize'],
    #             fontfamily=formatting['fontname'])
    
    # # Draw split
    # # proportion of vertical to horizontal extent of the slanted line
    # d = 0.5
    # kwargs = dict(marker=[(-d, -1), (d, 1)],
    #               linestyle='none',
    #               markersize=10,
    #               mec='k', mew=1,
    #               color='k', clip_on=False)
    # ax_left.plot([1], [0], transform=ax_left.transAxes, **kwargs)
    # ax_right.plot([0], [0], transform=ax_right.transAxes, **kwargs)
    
    # # Draw table
    # y_anchor = formatting['table_y_anchor'] * y_ticks[-1]
    # y_spacing = formatting['table_y_spacing'] * y_ticks[-1]
    # x_anchor = formatting['low_pCa_ticks'][0]
    
    # ax_right.text(x_anchor,
    #              y_anchor,
    #              'pCa$\\mathregular{_{50}}$',
    #              fontfamily=formatting['fontname'],
    #              fontsize=formatting['table_fontsize'],
    #              horizontalalignment='center',
    #              verticalalignment='center',
    #              clip_on=False)
    # ax_right.text(x_anchor - 1 * formatting['table_x_spacing'],
    #              y_anchor,
    #              'n$\\mathregular{_{H}}$',
    #              fontfamily=formatting['fontname'],
    #              fontsize=formatting['table_fontsize'],
    #              horizontalalignment='center',
    #              verticalalignment='center',
    #              clip_on=False)
    
    # # Add in data
    # for i in range(no_of_curves):
        
    #     pdf = pd_fit[i]
    
    #     # Deduce curve props        
    #     c_ind = i % len(formatting['marker_symbols'])
        
    #     y_anchor = y_anchor - y_spacing
    #     ax_right.plot(x_anchor + formatting['table_x_spacing'],
    #                 y_anchor,
    #                 formatting['marker_symbols'][c_ind],
    #                 mfc = formatting['color_set'][c_ind],
    #                 mec = formatting['color_set'][c_ind],
    #                 fillstyle = formatting['fill_styles'][c_ind],
    #                 markersize=formatting['marker_size'],
    #                 clip_on=False)
    
    #     ax_right.text(x_anchor,
    #                 y_anchor,
    #                 '%.2f' % pdf['pCa_50'],
    #                 fontfamily=formatting['fontname'],
    #                 fontsize=formatting['table_fontsize'],
    #                 horizontalalignment='center',
    #                 verticalalignment='center',
    #                 clip_on=False)
    #     ax_right.text(x_anchor - formatting['table_x_spacing'],
    #                 y_anchor,
    #                 '%.2f' % pdf['n_H'],
    #                 fontfamily=formatting['fontname'],
    #                 fontsize=formatting['table_fontsize'],
    #                 horizontalalignment='center',
    #                 verticalalignment='center',
    #                 clip_on=False)

    

def create_y_pCa_figure(fig_data, batch_file_string):
    """ Creates a y pCa figure based on dict fig_data """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    # Iterate through the folder structure looking for curve data

    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)

    # Set up current figure
    fig = plt.figure(constrained_layout=False)
    gs = fig.add_gridspec(nrows=1, ncols=6,
                          left=0.3, right=0.95, wspace=0.1,
                          bottom = 0.2)
    fig.set_size_inches([3.5, 2.5])
    ax_left = fig.add_subplot(gs[0,0])
    ax_right = fig.add_subplot(gs[0, 1:5])

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']

    # Store curve data in lists
    pCa_values = defaultdict(list)
    y_values = defaultdict(list)

    curve_counter = 1
    keep_going = True

    # Create lists to hold data
    curve_index = []
    hs_force = []
    hs_pCa = []
    hs_length = []
        
    # Keep track of max_y
    max_y = -np.inf
    
    # Set up a dictionary to store the curve data
    curve_data = dict()
    curve_data['curve'] = []

    while (keep_going):
        curve_folder = os.path.join(top_data_folder,
                                    ('%i' % curve_counter))

        if os.path.isdir(curve_folder):
            # Find the results files
            for file in os.listdir(curve_folder):
                if (file.endswith('.txt') and not file.startswith('rates')):
                    data_file_string = \
                        os.path.join(curve_folder, file)
                    d = pd.read_csv(data_file_string, delimiter='\t')
                    pCa_values[curve_counter-1].append(d['hs_1_pCa'].iloc[-1])
                    y = formatting['y_scaling_factor'] * \
                            d[fig_data['data_field']].iloc[-50:-1].mean() # take the mean force over last 50 points
                    y_values[curve_counter-1].append(y)
                    if (np.amax(y) > max_y):
                        max_y = np.amax(y)
                    
                    # Store data for subsequent output
                    curve_index.append(curve_counter)
                    hs_force.append(y)
                    hs_pCa.append(d['hs_1_pCa'].iloc[-1])
                    hs_length.append(d['hs_1_length'].iloc[-1])

            # Add in curve
            res=cv.fit_pCa_data(pCa_values[curve_counter-1],
                                y_values[curve_counter-1])
            
            # Store data
            d_parameters = pd.DataFrame({'pCa_50': res['pCa_50'],
                                        'n_H': res['n_H'],
                                        'y_min': res['y_min'],
                                        'y_amp': res['y_amp']},
                                        index=[0])
            d_fits = pd.DataFrame({'x_fit': res['x_fit'],
                                    'y_fit': res['y_fit']});
            
            d_curve = pd.concat([d_parameters, d_fits])
            curve_data['curve'].append(d_curve)
            
            # plot
            
            # Deduce curve props
            c_ind = (curve_counter-1) % len(formatting['marker_symbols'])

            # Normalized data_field if required
            if formatting['y_normalized_to_max']:
                res['y_fit'] = res['y_fit']/max(y_values[curve_counter-1])
                y_values[curve_counter-1] = [x/max(y_values[curve_counter-1]) for x in y_values[curve_counter-1]]                

            for a in [ax_left, ax_right]:
                a.plot(pCa_values[curve_counter-1],
                    y_values[curve_counter-1],
                    formatting['marker_symbols'][c_ind],
                    fillstyle=formatting['fill_styles'][c_ind],
                    markersize=formatting['marker_size'],
                    markerfacecolor = formatting['color_set'][c_ind],
                    markeredgewidth=formatting['marker_edge_width'],
                    markeredgecolor=formatting['color_set'][c_ind])

                if formatting['labels'] != []:
                    a.plot(res['x_fit'], res['y_fit'],'-',
                        color = formatting['color_set'][c_ind],
                        label = formatting['labels'][c_ind])               
                else:                   
                    a.plot(res['x_fit'], res['y_fit'],
                        color = formatting['color_set'][c_ind],
                        linestyle=formatting['line_styles'][c_ind])

            # Loop on to the next folder
            curve_counter = curve_counter + 1

        else:
            keep_going = False

    # Add legend if required

    if formatting['labels'] != []:
        ax_right.legend(loc=formatting['legend_location'],
                        fontsize = formatting['y_label_fontsize']-2)

    # Take lists and create a data frame
    r = pd.DataFrame({'curve': curve_index,
                      'hs_pCa': hs_pCa,
                      'hs_force': hs_force,
                      'hs_length': hs_length})

    # Tidy up axis
    ax_left.set_xlim(formatting['high_pCa_tick'] +
                     formatting['high_pCa_span']*np.asarray([0.5, -0.5]))
    right_lims = [formatting['low_pCa_ticks'][0] + formatting['low_pCa_pad'],
                  formatting['low_pCa_ticks'][-1]]
    ax_right.set_xlim(right_lims)

    ax_left.spines['right'].set_visible(False)
    ax_left.spines['top'].set_visible(False)
    ax_right.spines['left'].set_visible(False)
    ax_right.spines['right'].set_visible(False)
    ax_right.spines['top'].set_visible(False)
    ax_right.tick_params('y', left=False, labelleft=False)

    # Ticks
    ax_left.set_xticks([formatting['high_pCa_tick']])
    ax_left.set_ylabel(formatting['y_axis_label'],
                       loc='center',
                       verticalalignment='center',
                       labelpad=formatting['y_label_pad'],
                       fontfamily=formatting['fontname'],
                       fontsize=formatting['y_label_fontsize'],
                       rotation=formatting['y_label_rotation'])

    ax_right.set_xlabel('pCa',
                       loc='center',
                       verticalalignment='center',
                       labelpad=formatting['x_label_pad'],
                       fontfamily=formatting['fontname'],
                       fontsize=formatting['x_label_fontsize'])
    ax_right.set_xticks(formatting['low_pCa_ticks'])


    y_ticks = [0, ut.multiple_greater_than(max_y,
                                           0.05*np.power(10, np.ceil(np.log10(max_y))))]

    for a in [ax_left, ax_right]:
        a.set_xticklabels(a.get_xticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        a.set_ylim(y_ticks)
        a.set_yticks(y_ticks)
        a.set_yticklabels(a.get_yticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])

        if formatting['y_normalized_to_max']:
            y_ticks = [0,1]
            a.set_ylim([0,1.03])
            a.spines['left'].set_bounds(y_ticks)
            a.set_yticks(y_ticks)
            a.set_yticklabels(a.get_yticks(),
                fontsize=formatting['tick_labels_fontsize'],
                fontfamily=formatting['fontname'])

    # Draw split
    # proportion of vertical to horizontal extent of the slanted line
    d = 0.5
    kwargs = dict(marker=[(-d, -1), (d, 1)],
                  linestyle='none',
                  markersize=10,
                  mec='k', mew=1,
                  color='k', clip_on=False)
    ax_left.plot([1], [0], transform=ax_left.transAxes, **kwargs)
    ax_right.plot([0], [0], transform=ax_right.transAxes, **kwargs)

    # Draw table
    y_anchor = formatting['table_y_anchor'] * y_ticks[-1]
    y_spacing = formatting['table_y_spacing'] * y_ticks[-1]
    x_anchor = formatting['low_pCa_ticks'][0]

    ax_right.text(x_anchor,
                 y_anchor,
                 'pCa$\\mathregular{_{50}}$',
                 fontfamily=formatting['fontname'],
                 fontsize=formatting['table_fontsize'],
                 horizontalalignment='center',
                 verticalalignment='center',
                 clip_on=False)
    ax_right.text(x_anchor - 1 * formatting['table_x_spacing'],
                 y_anchor,
                 'n$\\mathregular{_{H}}$',
                 fontfamily=formatting['fontname'],
                 fontsize=formatting['table_fontsize'],
                 horizontalalignment='center',
                 verticalalignment='center',
                 clip_on=False)
    
    # Add in data
    for i,c in enumerate(curve_data['curve']):

        # Deduce curve props        
        c_ind = i % len(formatting['marker_symbols'])
        
        y_anchor = y_anchor - y_spacing
        ax_right.plot(x_anchor + formatting['table_x_spacing'],
                    y_anchor,
                    formatting['marker_symbols'][c_ind],
                    mfc = formatting['color_set'][c_ind],
                    mec = formatting['color_set'][c_ind],
                    fillstyle = formatting['fill_styles'][c_ind],
                    markersize=formatting['marker_size'],
                    clip_on=False)

        ax_right.text(x_anchor,
                    y_anchor,
                    '%.2f' % c['pCa_50'].iloc[0],
                    fontfamily=formatting['fontname'],
                    fontsize=formatting['table_fontsize'],
                    horizontalalignment='center',
                    verticalalignment='center',
                    clip_on=False)
        ax_right.text(x_anchor - formatting['table_x_spacing'],
                    y_anchor,
                    '%.2f' % c['n_H'].iloc[0],
                    fontfamily=formatting['fontname'],
                    fontsize=formatting['table_fontsize'],
                    horizontalalignment='center',
                    verticalalignment='center',
                    clip_on=False)

    # Save the figure
    if ('output_image_file' in fig_data):
        if (fig_data['relative_to'] == 'this_file'):
                fig_data['output_image_file'] = \
                    os.path.join(base_folder,
                                 fig_data['output_image_file'])

        # Check dir exists
        dir_name = os.path.dirname(fig_data['output_image_file'])
        if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
            
        for f in fig_data['output_image_formats']:
            ofs = '%s.%s' % (fig_data['output_image_file'], f)
            print('Saving pCa_figure to: %s' % ofs)
            fig.savefig(ofs, dpi=200, bbox_inches='tight')
            
    plt.close()

    # Save the data as an excel file if required in the batch file
    if('output_data_file_string' in fig_data):
        if (fig_data['relative_to'] == 'this_file'):
            output_file_string = os.path.join(base_folder,
                                              fig_data['output_data_file_string'])
        else:
            output_file_string = fig_data['output_data_file_string']

        print('Writing tension-pCa data to %s' % output_file_string)
        with pd.ExcelWriter(output_file_string, engine='openpyxl') as writer:
            r.to_excel(writer, sheet_name = 'simulation_data', index=False)
            for (i,c) in enumerate(curve_data['curve']):
                c.to_excel(writer,
                           index=False,
                           sheet_name = 'curve_%i' % (i+1))


def create_fv_and_power_figure(fig_data, batch_file_string):
    """ Creates an fv and power figure """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    # Find the data folders
    base_folder = os.path.dirname(batch_file_string)

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']


    # Find the length fit model 
    if (not 'length_fit_mode' in fig_data): # fit mode not specified, default is exponential
        length_fit_mode = 'exponential'
    else:
        length_fit_mode = fig_data['length_fit_mode']

    curve_counter = 1
    keep_going = True

    # Create lists to hold data
    curve = []
    release_index = []
    m_force = []
    m_force_isometric = []
    m_force_passive = []
    m_force_passive_corrected = []
    m_velocity = []
    m_velocity_l0_per_s = []
    m_power = []
    m_power_passive_corrected = []

    while keep_going:
        curve_folder = os.path.join(top_data_folder,
                                    ('%i' % curve_counter))
        if (os.path.isdir(curve_folder)):
            
            # Find the results files and sort them in natural order
            results_files = []
            file_ind = []
            for file in os.listdir(curve_folder):
                if file.endswith('.txt'):
                    results_files.append(file)
                    file_ind.append(int(file.split('_')[1]))
            si = np.argsort(np.asarray(file_ind))
            results_files = [results_files[i] for i in si]
            
            # Make a figure to check 
            fig = plt.figure(constrained_layout=True)
            gs = fig.add_gridspec(nrows=2, ncols=1,
                                  wspace = 0.5,
                                  hspace=0.1)
            fig.set_size_inches([3.5, 4])
            ax_for = fig.add_subplot(gs[0,0])
            ax_len = fig.add_subplot(gs[1, 0])
            
            for file in results_files:
                data_file_string = os.path.join(curve_folder, file)
                
                # Display, as this can be slow
                print('Fitting shortening velocity for: %s' %
                      data_file_string)

                # Load up the results file
                d = pd.read_csv(data_file_string, delimiter='\t')
                initial_ml = d['m_length'].iloc[0] # muscle length at t = 0

                # Filter to fit time_interval
                d_fit = d.loc[(d['time'] >= fig_data['fit_time_interval_s'][0]) &
                              (d['time'] <= fig_data['fit_time_interval_s'][-1])]
                
                # Filter for display
                d_display = d.loc[(d['time'] >= (fig_data['fit_time_interval_s'][0] -
                                                 (0.5 * (fig_data['fit_time_interval_s'][-1] -
                                                       fig_data['fit_time_interval_s'][0]))))]
                
                ax_for.plot(d_display['time'], d_display['m_force'], 'b-')
                ax_len.plot(d_display['time'], d_display['m_length'], 'b-')

                # Now do the fit
                if length_fit_mode == 'exponential':

                    # Set the time of clamp as t = 0 
                    if  ('sim_release_s' in  fig_data):

                        time_offset = d_fit['time'] - fig_data['sim_release_s']

                    # if time of clamp not specified, set start fitting time as t = 0

                    else: 
                        time_offset = d_fit['time'] - fig_data['fit_time_interval_s'][0]                            

                    # vel_data = cv.fit_exponential_decay(time_offset.to_numpy(),
                    #                         d_fit['hs_1_length'].to_numpy())
                    
                    vel_data = cv.fit_shortening_length_trace(time_offset.to_numpy(),
                                                              d_fit['m_length'].to_numpy())
                    
                    # Shortening velocity in ML s-1:
                    dydt = -(vel_data['b'] * vel_data['c'])
                    
                    m_vel = -1e-9 * dydt
                    m_vel_l0_per_s = -dydt / initial_ml
                        

                    # hs_vel = -1e-9* vel_data['amp'] * vel_data['k']  # velocity in m s^-1
                    # hs_vel_l0_per_s = -vel_data['amp'] * vel_data['k']/initial_hsl # velocity in ML s^-1

                    
                elif length_fit_mode == 'linear':

                    vel_data = cv.fit_straight_line(d_fit['time'].to_numpy(),
                                            d_fit['m_length'].to_numpy())

                    # Calculate velocity
                    m_vel = -1e-9*vel_data['slope'] # velocity in m s^-1
                    m_vel_l0_per_s = 1e9 * m_vel / initial_ml # velocity in ML s^-1
                   
                # Plot
                ax_len.plot(d_fit['time'], vel_data['y_fit'], 'r-')
          
                # Deduce isometric time
                d_isometric = d[(d['time'] > (fig_data['sim_release_s'] - 0.1)) &
                                (d['time'] < (fig_data['sim_release_s'] - 0.001))]
                m_for_isometric = d_isometric['m_force'].mean()
            
                # Get force and power

                m_for = d_fit['m_force'].mean()
                m_for_passive = d['m_force'][0]
                m_pow = m_vel * m_for / (1e-9 * initial_ml)
                m_pow_passive_corrected = \
                    m_vel * (m_for - m_for_passive) / \
                        (1e-9 * initial_ml)

                # Store data
                curve.append(curve_counter)
                release_index.append((os.path.basename(data_file_string)).split('_')[1])
                m_velocity.append(m_vel)
                m_velocity_l0_per_s.append(m_vel_l0_per_s)
                m_force.append(m_for)
                m_force_isometric.append(m_for_isometric)
                m_force_passive.append(m_for_passive)
                m_force_passive_corrected.append(
                    m_for - m_for_passive)
                m_power.append(m_pow)
                m_power_passive_corrected.append(m_pow_passive_corrected)
           
            
            fit_traces_string = ('fv_traces_%i' % curve_counter)
            if (fig_data['relative_to'] == 'this_file'):
                fit_traces_string = \
                        os.path.join(base_folder,
                                     fit_traces_string)
            else:
                dir_name = os.path.dirname(fig_data['output_image_file'])
                fit_traces_string = os.path.join(dir_name, fit_traces_string)

            # Check dir exists
            dir_name = os.path.dirname(fit_traces_string)
            
            if (not os.path.isdir(dir_name)):
                os.makedirs(dir_name)
                
            for f in fig_data['output_image_formats']:
                ofs = '%s.%s' % (fit_traces_string, f)
                print('Saving fit traces figure to: %s' % ofs)
                fig.savefig(ofs, dpi=200, bbox_inches='tight')

            curve_counter = curve_counter + 1

        else:
            keep_going = False

    # Take lists and create a data frame
    r = pd.DataFrame({'curve': curve,
                      'release_index': release_index,
                      'm_velocity': m_velocity,
                      'm_velocity_l0_per_s': m_velocity_l0_per_s,
                      'm_force': m_force,
                      'm_force_isometric': m_force_isometric,
                      'm_force_passive': m_force_passive,
                      'm_force_passive_corrected': m_force_passive_corrected,
                      'm_power': m_power,
                      'm_power_passive_corrected': m_power_passive_corrected})
    
    # Drop rows with NaNs, first replace empty with nan
    r.replace('', np.nan, inplace=True)
    r.dropna(inplace=True)
    
    # Now cycle through the curves fitting hybperbolas to each condition
    # This allows you to calculate force relative to isometric and
    # v relative to V_max
    
    for c in range(1, curve_counter):
        # First get the indices for the curve in the original dataframe        
        vi = r.index[r['curve'] == c].tolist()
        
        # Now make a new frame for just that curve
        rc = r[r['curve'] == c].copy()

        # Now calculate the relative force, relative velocity, and relative power            
        rc['m_f_to_f_max'] = rc['m_force'] / rc['m_force_isometric']
        rc['m_rel_power'] = rc['m_f_to_f_max'] * rc['m_velocity_l0_per_s']
        
        rc['m_f_to_f_max_pas_corrected'] = \
            (rc['m_force'] - rc['m_force_passive']) / \
                (rc['m_force_isometric'] - rc['m_force_passive'])
        
        # Now add these values back into the main frame
        r.loc[vi, 'm_f_to_f_max'] = rc['m_f_to_f_max']
        r.loc[vi, 'm_rel_power'] = rc['m_rel_power']
        r.loc[vi, 'm_f_to_f_max_pas_corrected'] = \
            rc['m_f_to_f_max_pas_corrected']
    

    # Create a figure
    if ('output_image_file' in fig_data):

        # Make a figure
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(nrows=2, ncols=2,
                              wspace = 0.1,
                              hspace=0.1)
        fig.set_size_inches([7, 5])
        ax_fv = fig.add_subplot(gs[0,0])
        ax_pow = fig.add_subplot(gs[1, 0])
        ax_rel_fv = fig.add_subplot(gs[0,1])
        ax_rel_pow = fig.add_subplot(gs[1,1])

        # Hold ticks
        f_ticks = np.asarray([])
        rel_f_ticks = np.asarray([])
        v_ticks = np.asarray([])
        rel_v_ticks = np.asarray([])
        p_ticks = np.asarray([])
        rel_p_ticks = np.asarray([])
        
        # Set up dict to store the curve data
        curve_data = dict()
        curve_data['curve'] = []
        
        # Cycle through curves
        for c in range(1, curve_counter):
            # Pull off the curve data
            rc = r[r['curve'] == c].copy()
            
            # Deduce curve props        
            c_ind = (c-1) % len(formatting['marker_symbols'])
            
            # Plot the force velocity curve
            ax_fv.plot(rc['m_force'], rc['m_velocity'],
                       formatting['marker_symbols'][c_ind],
                       formatting['marker_size'],
                       fillstyle=formatting['fill_styles'][c_ind],
                       color = formatting['color_set'][c_ind])
            
            #x,y  = cv.remove_outliers(rc['m_force'], rc['m_velocity'])
            #fv_curve = cv.fit_hyperbola(x, y)
            fv_curve = cv.fit_hyperbola(rc['m_force'], rc['m_velocity'])

            if formatting['labels'] != []:

                ax_fv.plot(fv_curve['x_fit'], fv_curve['y_fit'],
                           color=ax_fv.lines[-1].get_color(),
                           linestyle = formatting['line_styles'][c_ind],
                           label = formatting['labels'][c_ind])
            else:
                ax_fv.plot(fv_curve['x_fit'], fv_curve['y_fit'],
                    color=ax_fv.lines[-1].get_color(),
                    linestyle = formatting['line_styles'][c_ind],
                    label='Curve %i' % c)

            ax_pow.plot(rc['m_force'], rc['m_power'],
                       formatting['marker_symbols'][c_ind],
                       formatting['marker_size'],
                       fillstyle=formatting['fill_styles'][c_ind],
                        color = formatting['color_set'][c_ind])
            
            #x,y = cv.remove_outliers(rc['m_force'], rc['m_power'])
            #pow_curve = cv.fit_power_curve(x, y)
            pow_curve = cv.fit_power_curve(rc['m_force'], rc['m_power'])
            
            ax_pow.plot(pow_curve['x_fit'], pow_curve['y_fit'],
                        color=ax_pow.lines[-1].get_color(),
                        linestyle = formatting['line_styles'][c_ind])

            ax_rel_fv.plot(rc['m_f_to_f_max'], rc['m_velocity_l0_per_s'],
                           formatting['marker_symbols'][c_ind],
                           formatting['marker_size'],
                           fillstyle=formatting['fill_styles'][c_ind],
                           color = formatting['color_set'][c_ind])
            
            # x,y = cv.remove_outliers(rc['m_f_to_f_max'], rc['m_velocity_l0_per_s'])
            # rel_fv_curve = cv.fit_hyperbola(x, y)
            rel_fv_curve = cv.fit_hyperbola(rc['m_f_to_f_max'], rc['m_velocity_l0_per_s'])
            
            ax_rel_fv.plot(rel_fv_curve['x_fit'], rel_fv_curve['y_fit'],
                           color=ax_rel_fv.lines[-1].get_color(),
                           linestyle = formatting['line_styles'][c_ind])

            ax_rel_pow.plot(rc['m_f_to_f_max'], rc['m_rel_power'],
                            formatting['marker_symbols'][c_ind],
                            formatting['marker_size'],
                            fillstyle=formatting['fill_styles'][c_ind],
                            color = formatting['color_set'][c_ind])
            
            # x,y = cv.remove_outliers(rc['m_f_to_f_max'], rc['m_rel_power'])
            # rel_pow_curve = cv.fit_power_curve(x,y)
            rel_pow_curve = cv.fit_power_curve(rc['m_f_to_f_max'], rc['m_rel_power'])
            
            ax_rel_pow.plot(rel_pow_curve['x_fit'], rel_pow_curve['y_fit'],
                        color=ax_pow.lines[-1].get_color(),
                        linestyle = formatting['line_styles'][c_ind])

            # Store data to work out ticks later on
            f_ticks = np.concatenate((f_ticks, rc['m_force']))
            rel_f_ticks = np.concatenate((rel_f_ticks, rc['m_f_to_f_max']))
            v_ticks = np.concatenate((v_ticks, rc['m_velocity']))
            rel_v_ticks = np.concatenate((rel_v_ticks, rc['m_velocity_l0_per_s']))
            p_ticks = np.concatenate((p_ticks, rc['m_power']))
            rel_p_ticks = np.concatenate((rel_p_ticks, rc['m_rel_power']))
            
            # Deduce v_max
            fv_curve['v_max'] = ((fv_curve['x_0'] + fv_curve['a']) * 
                                    (fv_curve['b'] / fv_curve['a'])) - fv_curve['b']

            # Calculate v_max in muscle lengths per second
            fv_curve['v_max_l0_per_s'] = 1e9 * fv_curve['v_max'] / initial_ml

            
            # Store the curve data
            d_parameters = pd.DataFrame({'fv_x_0': fv_curve['x_0'],
                                         'fv_a': fv_curve['a'],
                                         'fv_b': fv_curve['b'],
                                         'fv_v_max': fv_curve['v_max'],
                                         'fv_v_max_l0_per_s': fv_curve['v_max_l0_per_s'],
                                         'pow_x_0': pow_curve['x_0'],
                                         'pow_a': pow_curve['a'],
                                         'pow_b': pow_curve['b'],
                                         'rel_fv_x_0': rel_fv_curve['x_0'],
                                         'rel_fv_a': rel_fv_curve['a'],
                                         'rel_fv_b': rel_fv_curve['b'],
                                         'rel_pow_x_0': rel_pow_curve['x_0'],
                                         'rel_pow_a': rel_pow_curve['a'],
                                         'rel_pow_b': rel_pow_curve['b']},
                                        index=[0])
            d_fits = pd.DataFrame({'fv_x_fit': fv_curve['x_fit'],
                                   'fv_y_fit': fv_curve['y_fit'],
                                   'pow_x_fit': pow_curve['x_fit'],
                                   'pow_y_fit': pow_curve['y_fit'],
                                   'rel_fv_x_fit': rel_fv_curve['x_fit'],
                                   'rel_fv_y_fit': rel_fv_curve['y_fit'],
                                   'rel_pow_x_fit': rel_pow_curve['x_fit'],
                                   'rel_pow_y_fit': rel_pow_curve['y_fit']})
            d_curve = pd.concat([d_parameters, d_fits])
            curve_data['curve'].append(d_curve)

        # Tidy up
            
        # Shortening velocity against stress
        xt = ut.tidy_limits(f_ticks)
        ax_fv.set_xlim(xt)
        ax_fv.set_xticks(xt)
        ax_fv.set_xticklabels(ax_fv.get_xticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        yt = ut.tidy_limits(v_ticks)
        ax_fv.set_ylim(yt)
        ax_fv.set_yticks(yt)
        ax_fv.set_yticklabels(ax_fv.get_yticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        ax_fv.set_xlabel('Force (N m$^{\\mathregular{-2}}$)',
                          fontfamily=formatting['fontname'],
                          loc='center')
        ax_fv.set_ylabel('Shortening\nvelocity\n(m s$^{\\mathregular{-1}}$)',
                          fontfamily=formatting['fontname'],
                          loc='center',
                          verticalalignment='center',
                          rotation=formatting['y_label_rotation'])
        ax_fv.legend()

        # Add legend on first panel if required

        if formatting['labels'] != []:
            ax_fv.legend(loc=formatting['legend_location'], fontsize = formatting['y_label_fontsize']-2)

        # Power in W m^-3 against Stress
        ax_pow.set_xlim(xt)
        ax_pow.set_xticks(xt)
        ax_pow.set_xticklabels(ax_pow.get_xticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        yt = ut.tidy_limits(p_ticks)
        ax_pow.set_ylim(yt)
        ax_pow.set_yticks(yt)
        ax_pow.set_yticklabels(ax_pow.get_yticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        ax_pow.set_xlabel('Force (N m$^{\\mathregular{-2}}$)',
                          fontfamily=formatting['fontname'],
                          loc='center')
        ax_pow.set_ylabel('Power\n(W m$^{\\mathregular{-3}}$)',
                          fontfamily=formatting['fontname'],
                          loc='center',
                          verticalalignment='center',
                          rotation=formatting['y_label_rotation'])
        
        # Rel velocity against rel force
        xt = ut.tidy_limits(rel_f_ticks)
        ax_rel_fv.set_xlim(xt)
        ax_rel_fv.set_xticks(xt)
        ax_rel_fv.set_xticklabels(ax_rel_fv.get_xticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        yt = ut.tidy_limits(rel_v_ticks)
        ax_rel_fv.set_ylim(yt)
        ax_rel_fv.set_yticks(yt)
        ax_rel_fv.set_yticklabels(ax_rel_fv.get_yticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        ax_rel_fv.set_xlabel('Relative force',
                          fontfamily=formatting['fontname'],
                          loc='center')
        ax_rel_fv.set_ylabel('Relative\nshortening\nvelocity\n(l$_{\\mathregular{0}}$ s$^{\\mathregular{-1}}$)',
                          fontfamily=formatting['fontname'],
                          loc='center',
                          verticalalignment='center',
                          rotation=formatting['y_label_rotation'])

        # Rel power against rel force
        xt = ut.tidy_limits(rel_f_ticks)
        ax_rel_pow.set_xlim(xt)
        ax_rel_pow.set_xticks(xt)
        ax_rel_pow.set_xticklabels(ax_rel_fv.get_xticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        yt = ut.tidy_limits(rel_p_ticks)
        ax_rel_pow.set_ylim(yt)
        ax_rel_pow.set_yticks(yt)
        ax_rel_pow.set_yticklabels(ax_rel_fv.get_yticks(),
                          fontsize=formatting['tick_labels_fontsize'],
                          fontfamily=formatting['fontname'])
        ax_rel_pow.set_xlabel('Relative force',
                          fontfamily=formatting['fontname'],
                          loc='center')
        ax_rel_pow.set_ylabel('Relative\npower',
                          fontfamily=formatting['fontname'],
                          loc='center',
                          verticalalignment='center',
                          rotation=formatting['y_label_rotation'])

    # Save the figure
    if ('output_image_file' in fig_data):
        if (fig_data['relative_to'] == 'this_file'):
                fig_data['output_image_file'] = \
                    os.path.join(base_folder,
                                 fig_data['output_image_file'])

        # Check dir exists
        dir_name = os.path.dirname(fig_data['output_image_file'])
        if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
            
        for f in fig_data['output_image_formats']:
            ofs = '%s.%s' % (fig_data['output_image_file'], f)
            print('Saving force_velocity_figure to: %s' % ofs)
            fig.savefig(ofs, dpi=200, bbox_inches='tight')
            
    plt.close()
        
    # Save the data as an excel file
    if ('output_data_file_string' in fig_data):
        if (fig_data['relative_to'] == 'this_file'):
            output_file_string = os.path.join(base_folder,
                                              fig_data['output_data_file_string'])
        else:
            output_file_string = fig_data['output_data_file_string']
    
        with pd.ExcelWriter(output_file_string, engine='openpyxl') as writer:
            r.to_excel(writer, sheet_name = 'simulation_data', index=False)
            if ('output_image_file' in fig_data):
                for (i,c) in enumerate(curve_data['curve']):
                    c.to_excel(writer,
                           index=False,
                           sheet_name = ('curve_%i' % (i+1)))


def create_ktr_figure(fig_data, batch_file_string):
    """ Creates ktr figure based on dict fig_data """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    # Iterate through the folder structure looking for curve data

    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']

    curve_counter = 1
    keep_going = True

    # Create lists to hold data
    curve = []
    hs_ktr = []
    hs_pCa = []
    hs_force = []

    while keep_going:
        curve_folder = os.path.join(top_data_folder,
                                    ('%i' % curve_counter))
        if (os.path.isdir(curve_folder)):
            # Find the results files
            for file in os.listdir(curve_folder):
                if file.endswith('.txt'):
                    data_file_string = os.path.join(curve_folder, file)

                    # Load up the results file
                    d = pd.read_csv(data_file_string, delimiter='\t')

                    # Filter to fit time_interval
                    d_fit = d.loc[(d['time'] >= fig_data['fit_time_interval_s'][0]) &
                                  (d['time'] <= fig_data['fit_time_interval_s'][-1])]

                    
                    # Set the time origin to 0
                    orig_time = d_fit['time'] - fig_data['fit_time_interval_s'][0]

                    ktr_data = cv.fit_exponential_recovery(orig_time.to_numpy(),
                                                    d_fit['force'].to_numpy())
                    
                    # Store data
                    curve.append(curve_counter)
                    hs_ktr.append(ktr_data['k'])
                    hs_pCa.append(d_fit['pCa'].iloc[-1])
                    hs_force.append(formatting['x_scaling_factor'] * d_fit['force'].iloc[-1])

            curve_counter = curve_counter + 1

        else:
            keep_going = False

    # Take lists and create a data frame
    r = pd.DataFrame({'curve': curve,
                    'hs_pCa': hs_pCa,
                    'hs_force': hs_force,
                    'hs_ktr': hs_ktr})

    # Create a figure
    if ('output_image_file_string' in fig_data):
        # Deduce the output file string
        if (fig_data['relative_to'] == 'this_file'):
            output_image_file_string = os.path.join(
                base_folder, fig_data['output_image_file_string'])
        else:
            output_image_file_string = fig_data['output_image_file_string']

        # Make a figure
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(nrows=1, ncols=1,
                              wspace = 0.5,
                              hspace=0.1)
        fig.set_size_inches([7, 3.5])
        ax_ktr = fig.add_subplot(gs[0,0])

        # Cycle through curves
        for c in range(1, curve_counter):
            rc = r[r['curve'] == c]

            if formatting['x_field']== "force": # plot ktr-force curve

                # Normalized data_field if required
                if formatting['x_normalized_to_max']:
                    rc['hs_force'] = rc['hs_force']/max(rc['hs_force'])

                # Plot the ktr curve

                if formatting['labels'] != []:
                    ax_ktr.plot(rc['hs_force'], rc['hs_ktr'], '--o', color = formatting['color_set'][c - 1], label = formatting['labels'][c - 1])
                else:
                    ax_ktr.plot(rc['hs_force'], rc['hs_ktr'], '--o', color = formatting['color_set'][c - 1])


                ax_ktr.set_xlabel(formatting['x_axis_label'],
                        fontfamily=formatting['fontname'],
                        loc='center')

                ax_ktr.set_ylabel('k$_{\\mathregular{tr}}$\n(s$^{\\mathregular{-1}}$)',
                          fontfamily=formatting['fontname'],
                          loc='center',
                          rotation=formatting['y_label_rotation'])


            if formatting['x_field'] == "pCa": # plot ktr-pCa curve 

                # Plot the ktr curve
                ax_ktr.plot(rc['hs_pCa'], rc['hs_ktr'], '--o', color = formatting['color_set'][c - 1])
                
                ax_ktr.set_xlabel('pCa',
                                    fontfamily=formatting['fontname'],
                                    loc='center')
                ax_ktr.set_ylabel('k$_{\\mathregular{tr}}$\n(s$^{\\mathregular{-1}}$)',
                                    fontfamily=formatting['fontname'],
                                    loc='center',
                                    rotation=formatting['y_label_rotation'])

                ax_ktr.invert_xaxis()

        # Add legend on ktr curve if required

        if formatting['labels'] != []:
            ax_ktr.legend(loc=formatting['legend_location'], fontsize = formatting['y_label_fontsize']-2)        
            
        # Save figure
        print('Saving ktr figure to: %s'% output_image_file_string)
        dir_name = os.path.dirname(output_image_file_string)
        if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
        fig.savefig(output_image_file_string)
        plt.close()

    # Save the data as an excel file
    if (fig_data['relative_to'] == 'this_file'):
        output_file_string = os.path.join(base_folder,
                                          fig_data['output_data_file_string'])
    else:
        output_file_string = fig_data['output_data_file_string']
    print('Writing k_tr data to %s' % output_file_string)
    r.to_excel(output_file_string,
               engine='openpyxl',
               index=False)

def superpose_ktr_plots(fig_data, batch_file_string):
    """ Superpose k_tr data from multiple result files """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']

    # Pull data
    # Loop through the txt files in top_data_folder
    # add to a list, then naturally sort the list

    results_files = []

    for file in os.listdir(top_data_folder): # iterate over all files in results_folder
        if file.endswith('.txt'): # this is a results file
            file = os.path.join(top_data_folder, file)
            results_files.append(file)

    results_files = natsorted(results_files) 

    # Create figure

    if ('output_image_file_string' in fig_data):
        # Deduce the output file string
        if (fig_data['relative_to'] == 'this_file'):
            output_image_file_string = os.path.join(
                base_folder, fig_data['output_image_file_string'])
        else:
            output_image_file_string = fig_data['output_image_file_string']

    fig = plt.figure(constrained_layout=False)
    fig.set_size_inches([3.5, 7])
    spec = gridspec.GridSpec(nrows=5, ncols=1, figure=fig,
                             wspace=1)
    axs=[]
    for r in range(0,5):
        axs.append(fig.add_subplot(spec[r,0]))

    # Plot data

    for i in range(0, len(results_files)):

        d = pd.read_csv(results_files[i], sep='\t')
        x = d['time']

        # Set max time on x axis

        max_time = ut.multiple_greater_than(max(d["time"]),
            0.05*np.power(10, np.ceil(np.log10(max(d["time"])))))

        # pCa

        axs[0].plot(x, d['pCa'], color = formatting['color_set'][i])

        # HSL length + command

        axs[1].plot(x, d['hs_length'], color = "black")
        axs[1].plot(x, d['hs_command_length'], color = "tab:orange", linestyle='--')

        if i == len(results_files) -1:
            axs[1].plot(x, d['hs_length'], color = "black", label = "Length")
            axs[1].plot(x, d['hs_command_length'], color = "tab:red", linestyle='--', label = "Length \ncommand")
            axs[1].legend(loc='upper left', bbox_to_anchor=[0.55, 0.9], fontsize = 12, handlelength=1.0)        

        # Force

        d['force'] = d['force']/1000

        if formatting['labels'] != []:        
            axs[2].plot(x, d['force'], label = formatting['labels'][i], color = formatting['color_set'][i], zorder=len(results_files) - i)

        else:
            axs[2].plot(x, d['force'], color = formatting['color_set'][i], zorder=len(results_files) - i)

        if formatting['labels'] != []:
            axs[2].legend(loc='upper left', bbox_to_anchor=[0.65, 0.65], fontsize = 11)

        # Actin pop
        axs[3].plot(x, d['a_pop_0'], '--' , color = formatting['color_set'][i], zorder=len(results_files) - i)
        axs[3].plot(x, d['a_pop_1'], '-' , color = formatting['color_set'][i], zorder=len(results_files) - i)

        # Myosin pop
        axs[4].plot(x, d['m_pop_0'], '--' , color = formatting['color_set'][i], zorder=len(results_files) - i)
        axs[4].plot(x, d['m_pop_1'], '-' , color = formatting['color_set'][i], zorder=len(results_files) - i)
        axs[4].plot(x, d['m_pop_2'], '-' , color = formatting['color_set'][i], zorder=len(results_files) - i)

    # Clean axis

    for i in range(5):
        
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['bottom'].set_visible(False)
        axs[i].set_xticks([])
        axs[i].set_xlim([0, max_time]) 
        
        for axis in ['top','bottom','left','right']:
            axs[i].spines[axis].set_linewidth(1.5)
            
        for tick in axs[i].yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
            
        axs[i].tick_params(direction = "out", length = 6, width = 1.5)

    axs[4].spines['bottom'].set_visible(True)

    axs[4].set_xlim([0, max_time]) 
    axs[4].set_xticks([0, max_time])
    axs[4].set_xlabel('Time (s)', labelpad = -5)
    
    # Y lables and axis limits

    axs[0].set_ylim([9, 4])
    axs[0].set_yticks([9, 4])
    axs[0].text(-0.2, 6.5, "pCa" , fontsize = 14, transform = axs[0].transData, ha='center', va='center')

    max_length =ut. multiple_greater_than(max(d["hs_command_length"]),
               0.001*np.power(10, np.ceil(np.log10(max(d["hs_command_length"])))))

    min_length = ut.multiple_greater_than(min(d["hs_command_length"]),
               0.001*np.power(10, np.ceil(np.log10(min(d["hs_command_length"])))))

    axs[1].set_ylim([min_length-0.5, max_length+0.5])
    axs[1].set_yticks([min_length, max_length])

    axs[1].text(-0.2, (max_length + min_length)/2, "HS length \n($\mathregular{\mu}$m)" , fontsize = 14, transform = axs[1].transData, ha='center', va='center')

    max_force = ut.multiple_greater_than(max(d["force"]),
               0.1*np.power(10, np.ceil(np.log10(max(d["force"])))))

    axs[2].set_ylim([0, max_force])
    axs[2].set_yticks([0, max_force])

    axs[2].text(-0.2, (0 + max_force)/2, 'Force \n(kN $\\mathregular{m}^{\mathregular{-2}}$)' , fontsize = 14, transform = axs[2].transData, ha='center', va='center')

    axs[3].text(-0.2, 0.5, 'Thin \n filament' , fontsize = 14, transform = axs[3].transData, ha='center', va='center')

    axs[3].set_ylim([0, 1])
    axs[3].set_yticks([0, 1])

    axs[4].text(-0.2, 0.5, 'Thick \n filament' , fontsize = 14, transform = axs[4].transData, ha='center', va='center')
    axs[4].set_ylim([0, 1])
    axs[4].set_yticks([0, 1])

    # Save figure
    print('Saving superposing plot figure to %s' % output_image_file_string)
    # Check folder exists and make it if not
    dir_name = os.path.dirname(os.path.abspath(
        output_image_file_string))
    if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
    fig.savefig(output_image_file_string, bbox_inches='tight')
    plt.close()

def dose_response(fig_data, batch_file_string):
    """ Plot data as a function of myotrope dose """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']

    # Get the myotrope dose list

    dose = fig_data["dose_list"]
    drug_effect = fig_data["drug_effect"] # define if it is an increasing/decreasing Hill curve
    dose_counter = 0
    keep_going = True

   # Create lists to hold data
    curve = []
    y_values = []
    IC_50 = []
    n_H = []

    # Keep track of max_y
    max_y = -np.inf

     # Make a figure

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(nrows=1, ncols=1)
    fig.set_size_inches([4, 3.5])
    ax = fig.add_subplot(gs[0,0])

    while keep_going:
        curve_folder = os.path.join(top_data_folder,
                                    ('%i' % dose_counter))
        if (os.path.isdir(curve_folder)):
            # Find the results files
            for file in os.listdir(curve_folder):
                if file.endswith('.txt'):

                    data_file_string = os.path.join(curve_folder, file)
                    d = pd.read_csv(data_file_string, delimiter='\t')

                    y = formatting['y_scaling_factor'] * \
                            d[fig_data['data_field']].iloc[-1]

                    if (np.amax(y) > max_y):
                        max_y = np.amax(y)

                    curve.append(dose_counter)
                    y_values.append(y)

            dose_counter = dose_counter + 1

        else:
            keep_going = False

    # Plot the dose response curve
    ax.plot(dose, y_values, formatting['marker_symbols'][0], markersize=formatting['marker_size'],
        markerfacecolor = formatting['color_set'][0],
        markeredgewidth=0.0)
    # Add fit
    res=cv.fit_IC_50(dose,y_values, drug_effect)

    # Store data
    IC_50.append(res['IC_50'])
    n_H.append(res['n_H']) 

    ax.plot(res['x_fit'], res['y_fit'], '-', color = formatting['color_set'][0])

    # Tidy up axis

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Log scale
    ax.set_xscale('log')

    # Ticks and labels

    ax.set_ylabel(formatting['y_axis_label'],
                       loc='center',
                       verticalalignment='center',
                       labelpad=formatting['y_label_pad'],
                       fontfamily=formatting['fontname'],
                       fontsize=formatting['y_label_fontsize'],
                       rotation=formatting['y_label_rotation'])

    ax.set_xlabel(formatting['x_axis_label'],
                       loc='center',
                       verticalalignment='center',
                       labelpad=formatting['x_label_pad'],
                       fontfamily=formatting['fontname'],
                       fontsize=formatting['x_label_fontsize'])


    y_ticks = [0, ut.multiple_greater_than(max_y,
                            0.1*np.power(10, np.ceil(np.log10(max_y))))]

    ax.set_ylim(y_ticks)
    ax.set_yticks(y_ticks)

    xticks = [0.01,0.1,1,10,100]
    ax.set_xticks(xticks)

    # Draw table
    y_anchor = formatting['table_y_anchor'] * y_ticks[-1]
    y_spacing = formatting['table_y_spacing'] * y_ticks[-1]
    x_anchor = xticks[-2]

    ax.text(x_anchor,
                 y_anchor,
                 'IC$\\mathregular{_{50}}$ ($\\mathregular{\mu}$M)',
                 fontfamily=formatting['fontname'],
                 fontsize=formatting['table_fontsize'],
                 horizontalalignment='center',
                 verticalalignment='center',
                 clip_on=False)
    ax.text(x_anchor + 150 * formatting['table_x_spacing'],
                 y_anchor,
                 'n$\\mathregular{_{H}}$',
                 fontfamily=formatting['fontname'],
                 fontsize=formatting['table_fontsize'],
                 horizontalalignment='center',
                 verticalalignment='center',
                 clip_on=False)
    
    # Add in data
    for i in range(len(IC_50)):
        y_anchor = y_anchor - y_spacing
        ax.text(x_anchor,
                    y_anchor,
                    '%.2f' % IC_50[i],
                    fontfamily=formatting['fontname'],
                    fontsize=formatting['table_fontsize'],
                    horizontalalignment='center',
                    verticalalignment='center',
                    clip_on=False)
        ax.text(x_anchor + 150 * formatting['table_x_spacing'],
                    y_anchor,
                    '%.2f' % n_H[i],
                    fontfamily=formatting['fontname'],
                    fontsize=formatting['table_fontsize'],
                    horizontalalignment='center',
                    verticalalignment='center',
                    clip_on=False)

    # Take lists and create a data frame
    r = pd.DataFrame({'curve': curve,
                    'dose': dose,
                    fig_data['data_field']: y_values})

    # Save figure
    if ('output_image_file_string' in fig_data):
        # Deduce the output file string
        if (fig_data['relative_to'] == 'this_file'):
            output_image_file_string = os.path.join(
                base_folder, fig_data['output_image_file_string'])
        else:
            output_image_file_string = fig_data['output_image_file_string']

    print('Saving dose response figure to: %s' % output_image_file_string)
    dir_name = os.path.dirname(output_image_file_string)
    if (not os.path.isdir(dir_name)):
        os.makedirs(dir_name)
    fig.savefig(output_image_file_string)
    plt.close()

    # Save the data as an excel file
    if (fig_data['relative_to'] == 'this_file'):
        output_file_string = os.path.join(base_folder,
                                          fig_data['output_data_file_string'])
    else:
        output_file_string = fig_data['output_data_file_string']
    print('Writing dose reponse data to %s' % output_file_string)
    r.to_excel(output_file_string,
               engine='openpyxl',
               index=False)
    
def create_rates_figure(fig_data, batch_file_string):
    """ Create a rates figure """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']

    # Store the data from the rate files
    rates_data = dict()
    rates_data['condition'] = []
    
    model_counter = 1
    keep_going = True
    
    # Set up for building arrays
    m_schemes = []
    c_schemes = []
    
    max_no_of_rates = 0

    # Loop through data folders
    while (keep_going):
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % model_counter))

        if os.path.isdir(condition_folder):
            for file in os.listdir(condition_folder):
                if (file == 'rates.json'):
                    
                    fs = os.path.join(condition_folder, file)
                    with open(fs, 'r') as f:
                        d = json.load(f, strict=False)
                    
                    # Scan through looking for myosins
                    for (i,m) in enumerate(d['FiberSim_rates']['myosin']):
                        # Read in the scheme
                        # Pull of some data
                        # Append the scheme to the holder as a Pandas
                        # dataframe
                        
                        df = pd.read_csv(StringIO(m['scheme']), sep='\t')
                        no_of_rates = len(df.columns) - 1
                        max_no_of_rates = np.amax([max_no_of_rates, no_of_rates])
                        m_schemes.append(df)
                    
                    # And now for mybpc
                    for (i,c) in enumerate(d['FiberSim_rates']['mybpc']):
                        # Same as for myosin above
                        df = pd.read_csv(StringIO(c['scheme']), sep='\t')
                        no_of_rates = len(df.columns) - 1
                        max_no_of_rates = max([max_no_of_rates, no_of_rates])
                        c_schemes.append(df)
                        
            # Increment the folder counter
            model_counter = model_counter + 1
        else:
            keep_going = False

    # Set-up the figure
    no_of_cols = 2
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = max_no_of_rates, ncols = no_of_cols,
                             wspace = 1, hspace = 0.5)
    fig.set_size_inches([2 * no_of_cols, 2.5 * no_of_rates])

    ax = []

    # Loop through both m and c schemes
    for i in range(2):
        if (i==0):
            s = m_schemes
        else:
            s = c_schemes

        # Create the color_map
        color_map = formatting['color_set']
        if (len(s) > len(color_map)):
            color_map = []
            for k,c in enumerate(np.linspace(0, 1, len(s))):
                color_map.append(cm.rainbow(k))
       
        # Loop through the rates making plots
        for j in range(max_no_of_rates):
            ax.append(fig.add_subplot(spec[j,i]))
       
        # Loop throught the isotypes and conditions
        for k, d in enumerate(s):
            # Loop through the rates
            for j in range(max_no_of_rates):
                
                plot_index = (i * max_no_of_rates) + j
                
                r_string = 'r_%i' % (j+1)
                if (r_string in d.columns):
                    ax[plot_index].plot(d['x'], np.log10(d[r_string]), '-',
                               color = color_map[k])
                    
                # Add title
                if (j==0):
                    if (i==0):
                        title_string = 'Myosin'
                    else:
                        title_string = 'MyBPC'
                    ax[plot_index].set_title(title_string)

                if (k == (len(s)-1)):
                    # Tidy up
                    x_ticks = np.asarray([d['x'].min(), 0, d['x'].max()])
                    ax[plot_index].set_xlim(x_ticks[[0, -1]])
                    ax[plot_index].set_ylabel('log10 %s' % r_string)
                    ax[plot_index].set_ylim([0, 3])
                    ax[plot_index].set_yticks(np.arange(0,4,1))

    if (fig_data['output_image_file']):
        if (fig_data['relative_to'] == 'this_file'):
            fig_data['output_image_file'] = \
                os.path.join(base_folder,
                             fig_data['output_image_file'])
        # Check dir exists
        dir_name = os.path.dirname(fig_data['output_image_file'])
        if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
        
        for f in fig_data['output_image_formats']:
            ofs = '%s.%s' % (fig_data['output_image_file'], f)
            print('Saving rates to %s' % ofs)
            fig.savefig(ofs, dpi=200, bbox_inches='tight')

    plt.close()

def create_superposed_traces_figure(fig_data, batch_file_string):
    """ Collates simulations and plots superpositions of each condition
        in a column """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    layout = default_layout()
    if ('layout' in fig_data):
        for entry in fig_data['layout']:
            layout[entry] = fig_data['layout'][entry]

    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']
        
    # Loop through data folders
    model_counter = 1
    keep_going = True

    while (keep_going):
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % model_counter))

        if os.path.isdir(condition_folder):
            model_counter = model_counter + 1
        else:
            keep_going = False

    # Hold the no_of_conditions
    no_of_conditions = model_counter-1

    # Keep track of max and mins
    min_hsl = np.inf
    max_hsl = -np.inf
    min_force = np.inf
    max_force = -np.inf
    min_sc_length = np.inf
    max_sc_length = -np.inf
    
    # Create the color_map
    color_map = formatting['color_set']
    if (no_of_conditions > len(color_map)):
        color_map = []
        for i,c in enumerate(np.linspace(0, 1, no_of_conditions)):
            color_map.append(cm.rainbow(c))

    # Create the figure
    for i in range(no_of_conditions):
        # Pull off the data files
        condition_folder = os.path.join(top_data_folder,
                                        ('%i' % (i+1)))
        file_counter = 1
        for file in os.listdir(condition_folder):
            if ((file.endswith('.txt')) and not file.endswith('rates.txt')):
                fs = os.path.join(condition_folder, file)
                # Force to numeric
                d = pd.read_csv(fs, sep='\t', na_values=['-nan(ind)'])
                
                # Deduce the number of half-sarcomeres
                pCa_names = [col for col in d if col.endswith('pCa')]
                no_of_half_sarcomeres = len(pCa_names)
                
                if ((i==0) and (file_counter == 1)):
                    # Check for series compliance to determine how to
                    # make the figure
                    print(d)
                    if (d['sc_extension'].max() > 0):
                        no_of_rows = 7
                        pCa_row = 1
                        hs_length_row = 2
                        sc_length_row = 3
                        force_row = 4
                        thin_filament_row = 5
                        thick_filament_row = 6
                        mybpc_row = 7
                    else:
                        no_of_rows = 6
                        pCa_row = 1
                        hs_length_row = 2
                        force_row = 3
                        thin_filament_row = 4
                        thick_filament_row = 5
                        mybpc_row = 6
                    
                    # Set-up the figure
                    fig = plt.figure(constrained_layout = False)
                    spec = gridspec.GridSpec(nrows=no_of_rows,
                                             ncols=no_of_conditions,
                                             figure=fig,
                                             wspace = layout['grid_wspace'],
                                             hspace = layout['grid_hspace'])
                    fig.set_size_inches([3 * no_of_conditions, 2 * no_of_rows])

                    ax=[]
                
                if ('x_ticks' in formatting):
                    d = d[(d['time'] > formatting['x_ticks'][0]) &
                          (d['time'] <= formatting['x_ticks'][-1])]

                # Keep track of max and mins
                for hs in range(no_of_half_sarcomeres):
                    hs_label = 'hs_%i_length' % (hs+1)
                    min_hsl = np.amin([min_hsl, d[hs_label].min()])
                    max_hsl = np.amax([max_hsl, d[hs_label].max()])
                    
                min_force = np.amin([min_force, d['hs_1_force'].min()])
                max_force = np.amax([max_force, d['hs_1_force'].max()])

                min_sc_length = np.amin([min_sc_length, d['sc_extension'].min()])
                max_sc_length = np.amax([max_sc_length, d['sc_extension'].max()])

                if (file_counter==1):
                    # Make the plots
                    for j in range(no_of_rows):
                        ax.append(fig.add_subplot(spec[j,i]))

                # Now plot
                plot_index = (i*no_of_rows) + (pCa_row - 1)
                ax[plot_index].plot(d['time'], d['hs_1_pCa'], '-',
                                    color = color_map[i],
                                    linewidth = formatting['data_linewidth'])
                if ('column_titles' in formatting):
                    ax[plot_index].title.set_text(formatting['column_titles'][i])
                    
                plot_index = (i*no_of_rows) + (hs_length_row - 1)
                ax[plot_index].plot(d['time'], d['m_length'], 'k-',
                                    linewidth = 1,
                                    label = "M length")
                ax[plot_index].plot(d['time'], d['hs_1_slack_length'], 'r-',
                                    linewidth = 1,
                                    label = "Slack length")
                for hs in range(no_of_half_sarcomeres):
                    hs_label = 'hs_%i_length' % (hs+1)
                    ax[plot_index].plot(d['time'], d[hs_label], '-',
                                    color = color_map[i],
                                    linewidth = formatting['data_linewidth'])
                # ax[plot_index].plot(d['time'], d['hs_1_command_length'], '--',
                #                     color = color_map[i],
                #                     linewidth = formatting['data_linewidth'])
                
                if (no_of_rows == 7):
                    plot_index = (i * no_of_rows) + (sc_length_row - 1)
                    ax[plot_index].plot(d['time'], d['sc_extension'], '-',
                                        color = color_map[i],
                                        linewidth = formatting['data_linewidth'])

                plot_index = (i * no_of_rows) + (force_row - 1)
                ax[plot_index].plot(d['time'], d['hs_1_force'], '-',
                                    color = color_map[i],
                                    linewidth = formatting['data_linewidth'],
                                    label='Total')
                
                for hs in range(no_of_half_sarcomeres):
                    titin_label = 'hs_%i_titin_force' % (hs+1)
                    plot_label = titin_label
                    if (hs>0):
                        plot_label = 'None'
                
                    ax[plot_index].plot(d['time'], d[titin_label], ':',
                                        color = color_map[i],
                                        linewidth = formatting['data_linewidth'],
                                        label=plot_label)

                    viscous_label = 'hs_%i_viscous_force' % (hs+1)
                    plot_label = viscous_label
                    if (hs>0):
                        plot_label = 'None'
                    
                    ax[plot_index].plot(d['time'], d[viscous_label], '--',
                                    color = color_map[i],
                                    linewidth = formatting['data_linewidth'],
                                    label=plot_label)

                # Now the a_states
                # Deduce the number of states
                a_pop_names = [col for col in d if ('_a_pop_' in col)]
                no_of_a_states = int(len(a_pop_names) / no_of_half_sarcomeres)
                
                plot_index = (i * no_of_rows) + (thin_filament_row - 1)
                for hs in range(no_of_half_sarcomeres):
                    for a_counter in range(no_of_a_states):
                        a_pop_string = ('hs_%i_a_pop_%i' %
                                        (hs+1, a_counter+1))
                        if ((file_counter == 1) and (hs==0)):
                            label = a_pop_string
                        else:
                            label = None
            
                        ax[plot_index].plot(d['time'], d[a_pop_string],
                                '-',
                                color = color_map[a_counter],
                                linewidth = formatting['data_linewidth'],
                                label=label)

                # Now the m_states
                # Deduce the number of states
                m_pop_names = [col for col in d if ('_m_pop_' in col)]
                no_of_m_states = int(len(m_pop_names) / no_of_half_sarcomeres)

                plot_index = (i * no_of_rows) + (thick_filament_row - 1)
                for hs in range(no_of_half_sarcomeres):
                    for m_counter in range(no_of_m_states):
                        m_pop_string = ('hs_%i_m_pop_%i' %
                                        (hs+1, m_counter+1))
                        if ((file_counter== 1) and (hs==0)):
                            label = m_pop_string
                        else:
                            label = None
                            
                        ax[plot_index].plot(d['time'], d[m_pop_string],
                                            '-',
                                            color = color_map[m_counter],
                                            linewidth = formatting['data_linewidth'],
                                            label=label)
                        
                # Now the c_states
                # Deduce the number of states
                c_pop_names = [col for col in d if ('_c_pop_' in col)]
                no_of_m_states = int(len(c_pop_names) / no_of_half_sarcomeres)

                plot_index = (i * no_of_rows) + (mybpc_row - 1)
                for hs in range(no_of_half_sarcomeres):
                    for c_counter in range(no_of_m_states):
                        c_pop_string = ('hs_%i_c_pop_%i' %
                                        (hs+1, c_counter+1))
                        if ((file_counter== 1) and (hs==0)):
                            label = c_pop_string
                        else:
                            label = None
                            
                        ax[plot_index].plot(d['time'], d[c_pop_string],
                                            '-',
                                            color = color_map[c_counter],
                                            linewidth = formatting['data_linewidth'],
                                            label=label)

                # Update counter
                file_counter = file_counter + 1

        # Handle formatting
        if (i==0):
            plot_index = (i * no_of_rows) + (pCa_row - 1)
            ax[plot_index].set_ylabel('pCa')

            plot_index = (i * no_of_rows) + (hs_length_row - 1)
            ax[plot_index].set_ylabel('HS length\n(nm)')
            
            if (no_of_rows == 7):
                plot_index = (i * no_of_rows) + (sc_length_row - 1)
                ax[plot_index].set_ylabel('Series\nextension\n(nm)')

            plot_index = (i * no_of_rows) + (force_row - 1)
            ax[plot_index].set_ylabel('Force\n(N m$^{-2}$)')

            plot_index = (i * no_of_rows) + (thin_filament_row - 1)
            ax[plot_index].set_ylabel('Thin\nfilament')

            plot_index = (i * no_of_rows) + (thick_filament_row - 1)
            ax[plot_index].set_ylabel('Thick\nfilament')

            plot_index = (i * no_of_rows) + (mybpc_row - 1)
            ax[plot_index].set_ylabel('MyBP-C')


    # Set limits
    for i in range(no_of_conditions):

        plot_index = (i * no_of_rows) + (pCa_row - 1)
        y_ticks = [4, 9]
        ax[plot_index].set_ylim(y_ticks)
        ax[plot_index].invert_yaxis()
        ax[plot_index].set_yticks(y_ticks)
        if ('superposed_x_ticks' in fig_data):
            ax[plot_index].set_xlim(fig_data['superposed_x_ticks'])
            ax[plot_index].set_xticks(fig_data['superposed_x_ticks'])

        plot_index = plot_index + 1
        y_ticks = [min_hsl-10, max_hsl+10]
        ax[plot_index].set_ylim(y_ticks)
        ax[plot_index].set_yticks(y_ticks)
        if ('superposed_x_ticks' in fig_data):
            ax[plot_index].set_xlim(fig_data['superposed_x_ticks'])
            ax[plot_index].set_xticks(fig_data['superposed_x_ticks'])
            
        if (no_of_rows == 7):
            plot_index = (i * no_of_rows) + (sc_length_row - 1)
            y_ticks = [min_sc_length-10, max_sc_length+10]
            ax[plot_index].set_ylim(y_ticks)
            ax[plot_index].set_yticks([y_ticks[0], 0, y_ticks[-1]])
            ax[plot_index].set_ylabel('Series\nextension\n(nm)')

        plot_index = (i * no_of_rows) + (force_row - 1)
        y_ticks = np.asarray([min_force, 0, max_force])
        ax[plot_index].set_ylim(y_ticks[[0, -1]])
        ax[plot_index].set_yticks(y_ticks)
        
        if ('superposed_x_ticks' in fig_data):
            ax[plot_index].set_xlim(fig_data['superposed_x_ticks'])
            ax[plot_index].set_xticks(fig_data['superposed_x_ticks'])

        plot_index = (i * no_of_rows) + (thin_filament_row - 1)
        y_ticks = [0, 1]
        ax[plot_index].set_ylim(y_ticks)
        ax[plot_index].set_yticks(y_ticks)
        ax[plot_index].legend(loc=formatting['legend_location'],
                      handlelength=formatting['legend_handlelength'],
                      bbox_to_anchor=(
                          formatting['legend_bbox_to_anchor'][0],
                          formatting['legend_bbox_to_anchor'][1]),
                      prop={'family': formatting['fontname'],
                        'size': formatting['legend_fontsize']})
        
        if ('superposed_x_ticks' in fig_data):
            ax[plot_index].set_xlim(fig_data['superposed_x_ticks'])
            ax[plot_index].set_xticks(fig_data['superposed_x_ticks'])


        plot_index = (i * no_of_rows) + (thick_filament_row - 1)
        y_ticks = [0, 1]
        ax[plot_index].set_ylim(y_ticks)
        ax[plot_index].set_yticks(y_ticks)
        ax[plot_index].legend(loc=formatting['legend_location'],
                      handlelength=formatting['legend_handlelength'],
                      bbox_to_anchor=(
                          formatting['legend_bbox_to_anchor'][0],
                          formatting['legend_bbox_to_anchor'][1]),
                      prop={'family': formatting['fontname'],
                        'size': formatting['legend_fontsize']})
        
        if ('superposed_x_ticks' in fig_data):
            ax[plot_index].set_xlim(fig_data['superposed_x_ticks'])
            ax[plot_index].set_xticks(fig_data['superposed_x_ticks'])

        plot_index = (i * no_of_rows) + (mybpc_row - 1)
        y_ticks = [0, 1]
        ax[plot_index].set_ylim(y_ticks)
        ax[plot_index].set_yticks(y_ticks)
        ax[plot_index].legend(loc=formatting['legend_location'],
                              handlelength=formatting['legend_handlelength'],
                              bbox_to_anchor=(
                                  formatting['legend_bbox_to_anchor'][0],
                                  formatting['legend_bbox_to_anchor'][1]),
                              prop={'family': formatting['fontname'],
                                'size': formatting['legend_fontsize']})
        
        if ('superposed_x_ticks' in fig_data):
            ax[plot_index].set_xlim(fig_data['superposed_x_ticks'])
            ax[plot_index].set_xticks(fig_data['superposed_x_ticks'])
            
    # Remove unnecessary spines
    for i in range(no_of_conditions * no_of_rows):
        ax[i].spines['top'].set_color('None')
        ax[i].spines['right'].set_color('None')

    if (fig_data['output_image_file']):
        if (fig_data['relative_to'] == 'this_file'):
            fig_data['output_image_file'] = \
                os.path.join(base_folder,
                             fig_data['output_image_file'])
        # Check dir exists
        dir_name = os.path.dirname(fig_data['output_image_file'])
        if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
        
        for f in fig_data['output_image_formats']:
            ofs = '%s.%s' % (fig_data['output_image_file'], f)
            print('Saving superposed traces to %s' % ofs)
            fig.savefig(ofs, dpi=200, bbox_inches='tight')

    plt.close()

def create_k_tr_analysis_figure(fig_data, batch_file_string):
    """ Creates a k_tr analysis figure """

    # Pull default formatting, then overwrite any values from
    # input file
    formatting = default_formatting()
    if ('formatting' in fig_data):
        for entry in fig_data['formatting']:
            formatting[entry] = fig_data['formatting'][entry]

    # Find the data folders
    base_folder = os.path.dirname(batch_file_string)

    if (fig_data['relative_to'] == 'this_file'):
        top_data_folder = os.path.join(base_folder,
                                       fig_data['results_folder'])
    else:
        top_data_folder = fig_data['results_folder']

    curve_counter = 1
    keep_going = True

    # Create lists to hold data
    curve = []
    force = []
    k_tr = []
    k_tr_amp = []
    pCa = []
    hs_length = []
    k_tr_r_squared = []

    # And dicts to hold traces
    sims = dict()
    sims['raw'] = []
    sims['fit'] = []

    while keep_going:
        curve_folder = os.path.join(top_data_folder,
                                    ('%i' % curve_counter))
        if (os.path.isdir(curve_folder)):
            # Find the results files
            for file in os.listdir(curve_folder):
                if (file.endswith('.txt') and not file.startswith('rates')):
                    data_file_string = os.path.join(curve_folder, file)

                    # Load up the results file
                    d = pd.read_csv(data_file_string, delimiter='\t')

                    # Filter to fit time_interval
                    d_fit = d.loc[(d['time'] >= fig_data['k_tr_fit_time_s'][0]) &
                                  (d['time'] <= fig_data['k_tr_fit_time_s'][-1])].copy()

                    # Pull off time offset
                    x = d_fit['time'].to_numpy()
                    x = x - x[0]
                    y = d_fit['hs_1_force'].to_numpy()
                    try:
                        k_tr_data = cv.fit_exponential_recovery(x,y)
                    except:
                        k_tr_data = dict()
                        k_tr_data['k'] = np.nan
                        k_tr_data['amp'] = np.nan
                        k_tr_data['x_fit'] = x
                        k_tr_data['y_fit'] = y

                    # Store some values
                    curve.append(curve_counter)
                    force.append(d['hs_1_force'].iloc[-1])
                    k_tr.append(k_tr_data['k'])
                    k_tr_amp.append(k_tr_data['amp'])
                    pCa.append(d['hs_1_pCa'].iloc[-1])
                    hs_length.append(d['hs_1_length'].iloc[-1])
                    k_tr_r_squared.append(k_tr_data['r_squared'])

                    d_fit['x_fit'] = k_tr_data['x_fit'] + d_fit['time'].iloc[0]
                    d_fit['y_fit'] = k_tr_data['y_fit']

                    sims['raw'].append(d)
                    sims['fit'].append(d_fit)

            curve_counter = curve_counter + 1

        else:
            keep_going = False

    # Make a dataframe from the lists
    r = pd.DataFrame({'curve': curve,
                      'pCa': pCa,
                      'force': force,
                      'hs_length': hs_length,
                      'k_tr': k_tr,
                      'k_tr_amp': k_tr_amp,
                      'k_tr_r_squared': k_tr_r_squared})

    # Save the data as an excel file if required in the batch file
    if(fig_data['output_data_file_string']):
        if (fig_data['relative_to'] == 'this_file'):
            output_file_string = os.path.join(base_folder,
                                              fig_data['output_data_file_string'])
        else:
            output_file_string = fig_data['output_data_file_string']

        print('Writing k_tr_analysis data to %s' % output_file_string)
        with pd.ExcelWriter(output_file_string, engine='openpyxl') as writer:
            r.to_excel(writer, sheet_name = 'simulation_data', index=False)
            # for (i,c) in enumerate(curve_data['curve']):
            #     c.to_excel(writer,
            #                 index=False,
            #                 sheet_name = 'curve_%i' % (i+1))

    # Create a figure
    if ('output_image_file' in fig_data):

        # Make a figure
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(nrows=2, ncols=2,
                              wspace = 0.5,
                              hspace=0.5)
        fig.set_size_inches([7, 5])
        ax_force = fig.add_subplot(gs[0,0])
        ax_hsl = fig.add_subplot(gs[1, 0])
        ax_k_tr_force = fig.add_subplot(gs[0,1])
        ax_k_tr_pCa = fig.add_subplot(gs[1,1])

        # Cycle through curves
        for (i,c) in enumerate(range(1, curve_counter)):
            vi = np.flatnonzero(r['curve'].to_numpy() == c)
            
            # Deduce curve props        
            c_ind = i % len(formatting['marker_symbols'])
            
            # Draw traces
            for v in vi:
                d_fit = sims['fit'][v]
                d_raw = sims['raw'][v]
                # Deduce drawing range
                t_draw_start_s = d_fit['time'].iloc[0] - \
                    (0.2 * (d_fit['time'].iloc[-1] - d_fit['time'].iloc[0]))
                d_raw = d_raw[d_raw['time'] > t_draw_start_s]

                ax_force.plot(d_raw['time'], d_raw['hs_1_force'],
                              color = formatting['color_set'][c_ind])
                ax_force.plot(d_fit['time'], d_fit['y_fit'], color='k')

                # Now the lengths
                ax_hsl.plot(d_raw['time'], d_raw['hs_1_command_length'],
                              color = 'k')
                ax_hsl.plot(d_raw['time'], d_raw['hs_1_length'],
                              color = formatting['color_set'][c_ind])

            # Now draw the plots
            r2 = r[r['curve'] == c]
            ax_k_tr_force.plot(r2['force'], r2['k_tr'],
                               formatting['marker_symbols'][c_ind],
                               markerfacecolor = formatting['color_set'][c_ind],
                               markeredgecolor = formatting['color_set'][c_ind],
                               fillstyle = formatting['fill_styles'][c_ind],
                               markeredgewidth=formatting['marker_edge_width'])

            if formatting['labels'] != []:
        
                ax_k_tr_force.plot(r2['force'], r2['k_tr'],
                    linestyle = formatting['line_styles'][c_ind],
                    color = formatting['color_set'][c_ind],
                    label = formatting['labels'][c_ind])
            else:
                ax_k_tr_force.plot(r2['force'], r2['k_tr'],
                    linestyle = formatting['line_styles'][c_ind],
                    color = formatting['color_set'][c_ind])

            ax_k_tr_pCa.plot(r2['pCa'], r2['k_tr'],
                               formatting['marker_symbols'][c_ind],
                               markerfacecolor = formatting['color_set'][c_ind],
                               markeredgecolor = formatting['color_set'][c_ind],
                               fillstyle = formatting['fill_styles'][c_ind],
                               markeredgewidth=formatting['marker_edge_width'])
            ax_k_tr_pCa.plot(r2['pCa'], r2['k_tr'],
                             linestyle = formatting['line_styles'][c_ind],
                             color = formatting['color_set'][c_ind])
        
        # Draw some labels
        ax_force.set_ylabel('Force (N m$^{-2}$)')
        ax_force.set_xlabel('Time (s)')

        ax_hsl.set_ylabel('Half-sarcomere length (nm)')
        ax_hsl.set_xlabel('Time (s)')

        ax_k_tr_force.set_xlabel('Force (N m$^{-2}$)')
        ax_k_tr_force.set_ylabel('k_tr (s$^{-1}$)')
        if ('k_tr_ticks' in fig_data):
            ax_k_tr_force.set_ylim(fig_data['k_tr_ticks'])

        ax_k_tr_pCa.set_xlabel('pCa')
        ax_k_tr_pCa.invert_xaxis()
        ax_k_tr_pCa.set_ylabel('k_tr (s$^{-1}$)')
        if ('k_tr_ticks' in fig_data):
            ax_k_tr_pCa.set_ylim(fig_data['k_tr_ticks'])
        
        # Add legend on ktr-force curve if required

        if formatting['labels'] != []:
            ax_k_tr_force.legend(loc=formatting['legend_location'], fontsize = formatting['y_label_fontsize']-2)

        # Save it
        if (fig_data['relative_to'] == 'this_file'):
                fig_data['output_image_file'] = \
                    os.path.join(base_folder,
                                 fig_data['output_image_file'])

        # Check dir exists
        dir_name = os.path.dirname(fig_data['output_image_file'])
        if (not os.path.isdir(dir_name)):
            os.makedirs(dir_name)
            
        for f in fig_data['output_image_formats']:
            ofs = '%s.%s' % (fig_data['output_image_file'], f)
            print('Saving k_tr_figure to: %s' % ofs)
            fig.savefig(ofs, dpi=200, bbox_inches='tight')

    plt.close()