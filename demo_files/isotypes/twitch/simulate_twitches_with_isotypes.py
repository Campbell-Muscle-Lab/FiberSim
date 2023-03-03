# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 12:20:56 2022

@author: kscamp3
"""

import os
import sys
import json
import copy

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path

def demo_twitches_with_isotypes():
    """ Characterize models """
    
    # Different isotypes
    rel_alpha = [0, 0.5, 1.0]

    # Variables
    base_setup_file_string = '../base/setup.json'
    base_model_file_string = '../base/base_model.json'

    batch_mode = 'characterize'
    
    generated_model_file_base = '../generated/model'
    generated_setup_file_string = '../generated/generated_setup.json'

    # Add protocols to the path
    sys.path.append('../../../code/fiberpy/fiberpy/package/modules/protocols')
    import protocols as prot

    # Generate a twitch protocol
    protocol_file_string = '../protocol/protocol_twitch.txt'
    
    tp = prot.create_twitch_protocol()
    prot.write_protocol_to_file(tp, protocol_file_string)
        
    # Load the base setup and base model
    with open(base_setup_file_string, 'r') as f:
        base_setup = json.load(f) 

    with open(base_model_file_string, 'r') as f:
        base_model = json.load(f)
        
    # Get the parent dir, because all paths have to be converted to absolutes
    base_dir = Path(base_setup_file_string).parent.absolute()
        
    # Loop through mod factors
    rep_model_file_strings = list()
    
    for i in range(len(rel_alpha)):

        # Copy the base model
        rep_model = copy.deepcopy(base_model)

#         # prop_dead_heads
#         rep_model['m_parameters']['m_isotype_proportions'] = \
#             [1.0-rel_dead_heads[i], rel_dead_heads[i]]
               
#         # force-dependence
#         y = np.asarray(rep_model['m_kinetics'][0]['scheme'][0]['transition'][0]['rate_parameters'],
#                  dtype=np.float32)
#         y[1] = y[1] * rel_force_dependence[i]
#         rep_model['m_kinetics'][0]['scheme'][0]['transition'][0]['rate_parameters'] = y.tolist()
            
#         # Generate a model file name
#         rfs = ('%s_%i.json' % (generated_model_file_base, i+1))
        
#         # Correct path
#         rfs = str(Path(os.path.join(base_dir, rfs)).resolve())
        
#         # Check the path exists
#         parent_dir = Path(rfs).parent.absolute()
#         if not os.path.isdir(parent_dir):
#             os.makedirs(parent_dir)

#         # Write the model
#         with open(rfs, 'w') as f:
#             json.dump(rep_model, f, indent=4)
            
#         # Convert to absolute path
#         rfs = str(Path(os.path.join(base_dir, rfs)).resolve())

#         # Add in to array
#         rep_model_file_strings.append(rfs)
    
#     # Now copy the setup
#     generated_setup = copy.deepcopy(base_setup)
    
#     # And change the file names to absolute paths
#     generated_setup['FiberSim_characterization']['model']['relative_to'] = 'false'
#     generated_setup['FiberSim_characterization']['model']['model_files'] = rep_model_file_strings
#     generated_setup['FiberSim_characterization']['model']['options_file'] = \
#         str(Path(os.path.join(base_dir,
#                               generated_setup['FiberSim_characterization']['model']['options_file'])).resolve())
    
#     # Loop through the characterizations, changing these paths
#     characterize = generated_setup['FiberSim_characterization']['characterization']
#     for i,c in enumerate(characterize):
#         c['relative_to'] = 'false'
#         c['sim_folder'] = str(Path(os.path.join(base_dir,
#                                                 c['sim_folder'])).resolve())
#         generated_setup['FiberSim_characterization']['characterization'][i] = c

#     # And finally the setup file
#     generated_setup_file_string = str(Path(os.path.join(base_dir,
#                                                         generated_setup_file_string)).resolve())
    
#     print("\n\n%s\n\n" % generated_setup_file_string)
    
#     with open(generated_setup_file_string, 'w') as f:
#         json.dump(generated_setup, f, indent=4)
        
#     # Generate a command line
#     cs = 'pushd \"%s\" & python FiberPy.py %s %s & popd' % \
#             (FiberSim_code_dir, batch_mode, generated_setup_file_string)
    
#     # And run it
#     os.system(cs)
    
# def create_figure():
    
#     color_set = ['k','b','r']
#     n_models = 3
#     n_lengths = 3
    
#     pCa_values = [8.0, 4.5]
#     n_pCa = len(pCa_values)
    
#     max_hs_force = 50000
    
#     marker_size = 4
#     line_width = 0.5
    
#     label_pad = 40
    
#     # Create figure
#     fig = plt.figure(constrained_layout=False)
#     fig.set_size_inches([9, 9])
#     spec = gridspec.GridSpec(nrows=5,
#                              ncols=3,
#                              figure=fig,
#                              wspace=2,
#                              hspace=1)
    
#     isometric_data_file_string = '../sim_data/pCa_length_control/sim_output/pCa_analysis.xlsx'
#     d = pd.read_excel(isometric_data_file_string)
    
#     ax=[]
    
#     for i in range(n_models):
#         for j in range(n_lengths):
#             ind = (i * n_lengths) + j + 1
#             d2 = d[d['curve'] == ind]
        
#             x = d2['hs_pCa']
#             y = d2['hs_force']
        
#             # Fit curve
#             res=cv.fit_pCa_data(x, y)
        
#             # Store data
#             d_parameters = pd.DataFrame({'pCa_50': res['pCa_50'],
#                                     'n_H': res['n_H'],
#                                     'y_min': res['y_min'],
#                                     'y_amp': res['y_amp']},
#                                     index=[0])
#             d_fits = pd.DataFrame({'x_fit': res['x_fit'],
#                                    'y_fit': res['y_fit']});

#             if (j==0):
#                 ax.append(fig.add_subplot(spec[0,i]))
                
#             ax[i].plot(x, y, 'o', color = color_set[i],
#                          markerfacecolor='none',
#                          markersize=marker_size)
#             ax[i].plot(res['x_fit'], res['y_fit'], '-',
#                          color=color_set[i],
#                          linewidth=line_width)
        
#         ax[i].invert_xaxis()
#         ax[i].set_ylim([0, max_hs_force])
#         ax[i].set_xlabel('pCa')
    
#     ax[0].set_ylabel('Stress\n(N m^-2)',
#                      rotation=0,
#                      labelpad = label_pad)

#     hsl =np.NaN * np.ones([n_pCa, n_models, n_lengths])
#     abs_f = np.NaN * np.ones([n_pCa, n_models, n_lengths])
#     rel_abs_f = np.NaN * np.ones([n_pCa, n_models, n_lengths])
#     act_f = np.NaN * np.ones([n_models, n_lengths])
#     rel_act_f = np.NaN * np.ones([n_models, n_lengths])
   
#     for (pCa_counter, pCa) in enumerate(pCa_values):
        
#         for i in range(n_models):
            
#             for j in range(n_lengths):
#                 ind = (i * n_lengths) + j + 1
#                 d2 = d[(d['curve'] == ind) & (d['hs_pCa'] == pCa)]
                
#                 hsl[pCa_counter, i, j] = d2['hs_length'].to_numpy()   
#                 abs_f[pCa_counter, i, j] = d2['hs_force'].to_numpy()
#                 rel_abs_f[pCa_counter, i, j] = abs_f[pCa_counter, i, j] / \
#                                             abs_f[pCa_counter, i, 0]
                                            
#                 if (pCa_counter==1):
#                     act_f[i,j] = abs_f[1,i,j] - abs_f[0,i,j]
#                     rel_act_f[i,j] = act_f[i,j] / act_f[i,0]
                
#             # Fit straight line
#             abs_fit = cv.fit_straight_line(hsl[pCa_counter,i,:], abs_f[pCa_counter,i,:])
#             abs_x_fit = np.linspace(hsl[pCa_counter,i,0], hsl[pCa_counter,i,-1], 50)
#             abs_y_fit = abs_fit['intercept'] + abs_fit['slope'] * abs_x_fit
            
#             rel_abs_fit = cv.fit_straight_line(hsl[pCa_counter,i,:], rel_abs_f[pCa_counter,i,:])
#             rel_abs_x_fit = np.linspace(hsl[pCa_counter,i,0], hsl[pCa_counter,i,-1], 50)
#             rel_abs_y_fit = rel_abs_fit['intercept'] + rel_abs_fit['slope'] * rel_abs_x_fit

#             if (i==0):
#                 ax.append(fig.add_subplot(spec[pCa_counter+1,0]))
#                 ax_abs = len(ax)-1
                
#                 ax.append(fig.add_subplot(spec[pCa_counter+1,1]))
#                 ax_rel = len(ax)-1

#             ax[ax_abs].plot(hsl[pCa_counter,i,:], abs_f[pCa_counter,i,:], 'o', color=color_set[i],
#                                         markerfacecolor='none',
#                                         markersize=marker_size)
#             ax[ax_abs].plot(abs_x_fit, abs_y_fit, '-',
#                                         color = color_set[i],
#                                         linewidth=line_width)
            
#             ax[ax_rel].plot(hsl[pCa_counter,i,:], rel_abs_f[pCa_counter,i,:], 'o',
#                                         color=color_set[i],
#                                         markerfacecolor='none',
#                                         markersize=marker_size)
#             ax[ax_rel].plot(rel_abs_x_fit, rel_abs_y_fit, '-',
#                                         color = color_set[i],
#                                         linewidth = line_width)
            
#             ax[ax_abs].set_ylim([0, max_hs_force])
#             ax[ax_abs].set_ylabel('pCa %.1f stress\n(N m^-2)' % pCa,
#                                   rotation=0,
#                                   labelpad = label_pad)
#             ax[ax_abs].set_xlabel('Half-sarcomere length (nm)')

#             ax[ax_rel].set_ylim([0, np.amax(rel_abs_y_fit)])
#             ax[ax_rel].set_ylabel('Relative\npCa %.1f stress\n(N m^-2)' % pCa,
#                                   rotation=0,
#                                   labelpad = label_pad)
#             ax[ax_rel].set_xlabel('Half-sarcomere length (nm)')

            
#             if (pCa_counter==1):
                
#                 act_fit = cv.fit_straight_line(hsl[pCa_counter,i,:], act_f[i,:])
#                 act_x_fit = np.linspace(hsl[pCa_counter,i,0], hsl[pCa_counter,i,-1], 50)
#                 act_y_fit = act_fit['intercept'] + act_fit['slope'] * act_x_fit
                
#                 rel_act_fit = cv.fit_straight_line(hsl[pCa_counter,i,:], rel_act_f[i,:])
#                 rel_act_x_fit = np.linspace(hsl[pCa_counter,i,0], hsl[pCa_counter,i,-1], 50)
#                 rel_act_y_fit = rel_act_fit['intercept'] + rel_act_fit['slope'] * rel_act_x_fit

#                 if (i==0):
#                     ax.append(fig.add_subplot(spec[pCa_counter+2,0]))
#                     ax_abs_act = len(ax)-1
                    
#                     ax.append(fig.add_subplot(spec[pCa_counter+2,1]))
#                     ax_rel_act = len(ax)-1

#                 ax[ax_abs_act].plot(hsl[pCa_counter,i,:], act_f[i,:], 'o',
#                                           color=color_set[i],
#                                           markerfacecolor='none',
#                                           markersize=marker_size)
#                 ax[ax_abs_act].plot(act_x_fit, act_y_fit, '-',
#                                         color = color_set[i],
#                                         linewidth = line_width)


#                 ax[ax_rel_act].plot(hsl[pCa_counter,i,:], rel_act_f[i,:], 'o',
#                                         color=color_set[i],
#                                         markerfacecolor='none',
#                                         markersize=marker_size)
#                 ax[ax_rel_act].plot(rel_act_x_fit, rel_act_y_fit, '-',
#                                         color = color_set[i],
#                                         linewidth = line_width)
                
#                 ax[ax_abs_act].set_ylim([0, max_hs_force])
#                 ax[ax_abs_act].set_ylabel('Active stress\n(N m^-2)',
#                                       rotation=0,
#                                       labelpad=label_pad)
#                 ax[ax_abs_act].set_xlabel('Half-sarcomere length (nm)')
                
#                 ax[ax_rel_act].set_ylim([0, np.amax(rel_abs_y_fit)])
#                 ax[ax_rel_act].set_ylabel('Relative\nactive stress\n(N m^-2)',
#                                       rotation=0,
#                                       labelpad = label_pad)
#                 ax[ax_rel_act].set_xlabel('Half-sarcomere length (nm)')


#     # Power
#     power_data_file_string = '../sim_data/fv_pCa45/isotonic/sim_output/fv_analysis.xlsx'
#     d = pd.read_excel(power_data_file_string)
    
#     for i in range(n_models):
        
#         d2 = d[d['curve'] == (i+1)]
#         abs_f = d2['hs_force'].to_numpy()
#         rel_f = abs_f / np.amax(abs_f)
#         p = d2['hs_power'].to_numpy()
        
#         abs_p_fit = cv.fit_power_curve(abs_f, p)
        
#         rel_p_x = abs_p_fit['x_fit'] / abs_p_fit['x_0']
        
#         if (i==0):
#             ax.append(fig.add_subplot(spec[4,0]))
#             ax_abs_pow = len(ax)-1
            
#             ax.append(fig.add_subplot(spec[4,1]))
#             ax_rel_pow = len(ax)-1
            
        
        
#         ax[ax_abs_pow].plot(abs_f, p, 'o', color=color_set[i],
#                       markerfacecolor='none',
#                       markersize=marker_size)
#         ax[ax_abs_pow].plot(abs_p_fit['x_fit'], abs_p_fit['y_fit'], '-', color=color_set[i],
#                       linewidth=line_width)
            
#         ax[ax_rel_pow].plot(rel_f, p, 'o', color=color_set[i],
#                       markerfacecolor='none',
#                       markersize=marker_size)
#         ax[ax_rel_pow].plot(rel_p_x, abs_p_fit['y_fit'], '-', color=color_set[i],
#                       linewidth=line_width)
        
#         ax[ax_abs_pow].set_ylabel('Power\n(W m^-3)',
#                               rotation=0,
#                               labelpad = label_pad)
#         ax[ax_abs_pow].set_xlabel('Stress\n(kN m^-2)')
        
#         ax[ax_rel_pow].set_ylabel('Power\n(W m^-3)',
#                               rotation=0,
#                               labelpad = label_pad)
#         ax[ax_rel_pow].set_xlabel('Relative Stress')

    
#     fig.align_labels()

#     output_file_string = '../sim_data/pCa_length_control/sim_output/summary_figure.png'
#     fig.savefig(output_file_string, dpi=200, bbox_inches='tight')
    

if __name__ == "__main__":
    demo_twitches_with_isotypes()