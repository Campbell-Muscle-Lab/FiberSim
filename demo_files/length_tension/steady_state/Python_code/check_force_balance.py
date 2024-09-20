import os
import sys
import json

from pathlib import Path

import numpy as np

import re

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



# Add the validation package

this_path = os.path.dirname(__file__)
code_path = os.path.join(this_path, '../../../../code/FiberPy/fiberpy/package')

code_path = str(Path(code_path).absolute().resolve())

sys.path.append(code_path)

from modules.half_sarcomere import half_sarcomere

def check_force_balance(dump_file_string, analysis_counter):
    """ Checks force balance on the dump file """
    
    hs = half_sarcomere.half_sarcomere(dump_file_string)
    
    # Work out size of problem
    no_of_thin_filaments = len(hs['thin'])
    no_of_thick_filaments = len(hs['thick'])
    
    a_nodes_per_thin_filament = hs['hs_data']['a_nodes_per_thin_filament']
    m_nodes_per_thick_filament = hs['hs_data']['m_nodes_per_thick_filament']
    
    hs_total_nodes = (no_of_thin_filaments * a_nodes_per_thin_filament) + \
        (no_of_thick_filaments * m_nodes_per_thick_filament)
    
    a_k_stiff = hs['thin'][0]['a_k_stiff']
    a_rl = hs['thin'][0]['a_inter_bs_rest_length']
    a_bs_per_node = hs['thin'][0]['a_bs_per_node']
    
    m_k_stiff = hs['thick'][0]['m_k_stiff']
    m_rl = hs['thick'][0]['m_inter_crown_rest_length']
    m_cbs_per_node = hs['thick'][0]['m_cbs_per_node']
    
    force_diff = np.zeros(hs_total_nodes)
    x_dev = np.zeros(hs_total_nodes)
    
    thin_ext = np.zeros((a_nodes_per_thin_filament-1, no_of_thin_filaments))
    thick_ext = np.zeros((m_nodes_per_thick_filament-1, no_of_thick_filaments))
    
    # Loop through thin filaments
    counter = 0
    for a_f, thin in enumerate(hs['thin']):
        x_node = thin['bs_x'][0::a_bs_per_node]
        
        for (i,x) in enumerate(x_node):
            # Pull off positions
            x = x_node[i]
            link_multiplier = 1
            if (i==0):
                x_neg = 0
                x_pos = x_node[i+1]
            elif (i == (len(x_node) - 1)):
                x_neg = x_node[i-1]
                x_pos = x
                link_multiplier = 0
            else:
                x_neg = x_node[i-1]
                x_pos = x_node[i+1]

            fd = link_multiplier * (a_k_stiff * (x_pos - x - a_rl)) - \
                    (a_k_stiff * (x - x_neg - a_rl))
        
            force_diff[counter] = fd
            
            counter = counter + 1
            
    for m_f, thick in enumerate(hs['thick']):
        x_node = thick['cb_x'][0::m_cbs_per_node]
        
        for (i,x) in enumerate(x_node):
            # Pull off positions
            x = x_node[i]
            link_multiplier = 1
            if (i==0):
                x_pos = hs['hs_data']['hs_length'] - thick['m_lambda']
                x_neg = x_node[i+1]
            elif (i == (len(x_node) - 1)):
                x_pos = x_node[i-1]
                x_neg = x
                link_multiplier = 0
            else:
                x_pos = x_node[i-1]
                x_neg = x_node[i+1]
        
            fd = (m_k_stiff * (x_pos - x - m_rl)) - \
                    (link_multiplier * m_k_stiff * (x - x_neg - m_rl))
                    
            force_diff[counter] = fd
            
            counter = counter + 1
 
    # Now add in corrections for titin
    for m_f, thick in enumerate(hs['thick']):
        
        # Get the position on the thick filament
        m_n = (hs['titin']['t_attach_m_node'] - 1) * m_cbs_per_node
        x_m = hs['thick'][m_f]['cb_x'][m_n]
        
        # Index for the m_row
        m_row = (no_of_thin_filaments * a_nodes_per_thin_filament) + \
            (m_f * m_nodes_per_thick_filament) + \
            (hs['titin']['t_attach_m_node'] - 1)
        
        for (i, thin_id) in enumerate(thick['nearest_actin_filaments']):
            a_n = (hs['titin']['t_attach_a_node'] - 1) * a_bs_per_node
            x_a = hs['thin'][thin_id]['bs_x'][a_n]
           
            fd = (hs['titin']['t_k_stiff'] * \
                      (x_m - x_a - hs['titin']['t_offset']))
                
            if (hs['titin']['t_passive_mode'] == 'exponential'):
                fd = fd + \
                    (np.sign(x_m - x_a) * hs['titin']['t_sigma'] * \
                     (np.exp(np.fabs(x_m - x_a - hs['titin']['t_offset']) /
                             hs['titin']['t_L']) - 1))
            
            a_row = (thin_id * a_nodes_per_thin_filament) + \
                (hs['titin']['t_attach_a_node'] -1)
            
            force_diff[a_row] = force_diff[a_row] + fd
            force_diff[m_row] = force_diff[m_row] - fd
                
    # Now add in correction for myosin
    for m_f, thick in enumerate(hs['thick']):
        for m_n, bound_to_a_f in enumerate(thick['cb_bound_to_a_f']):
            if (bound_to_a_f >= 0):
                a_f = bound_to_a_f
                a_n = thick['cb_bound_to_a_n'][m_n]
                
                x_m = thick['cb_x'][m_n]
                x_a = hs['thin'][a_f]['bs_x'][a_n]
                
                a_row = (a_f * a_nodes_per_thin_filament) + \
                    int(np.floor(a_n / thin['a_bs_per_node']))
                    
                m_row = (no_of_thin_filaments * a_nodes_per_thin_filament) + \
                    (m_f * m_nodes_per_thick_filament) + \
                    int(np.floor(m_n / thick['m_cbs_per_node']))
                
                ext = hs['hs_data']['cb_extensions'][thick['cb_state'][m_n]-1]
                
                fd = thick['m_k_cb'] * (x_m - x_a + ext)
                
                force_diff[a_row] = force_diff[a_row] + fd
                force_diff[m_row] = force_diff[m_row] - fd
                
    # And finally correction for mybpc
    for m_f, thick in enumerate(hs['thick']):
        for m_n, bound_to_a_f in enumerate(thick['pc_bound_to_a_f']):
            if (bound_to_a_f >= 0):
                a_f = bound_to_a_f
                a_n = thick['pc_bound_to_a_n'][m_n]
                
                thick_node_index = thick['pc_node_index'][m_n]
                x_m = thick['cb_x'][thick_node_index * thick['m_cbs_per_node']]
                x_a = hs['thin'][a_f]['bs_x'][a_n]
                
                a_row = (a_f * a_nodes_per_thin_filament) + \
                         int(np.floor(a_n / thin['a_bs_per_node']))
                
                m_row = (no_of_thin_filaments * a_nodes_per_thin_filament) + \
                    (m_f * m_nodes_per_thick_filament) + \
                        thick_node_index
                        
                fd = 0 * thick['c_k_stiff'] * (x_m - x_a)
                
                force_diff[a_row] = force_diff[a_row] + fd
                force_diff[m_row] = force_diff[m_row] - fd
        
    # Calculate deviations
    for i in range(hs_total_nodes):
        if (i < (no_of_thin_filaments * a_nodes_per_thin_filament)):
            x_dev[i] = force_diff[i] / a_k_stiff
        else:
            x_dev[i] = force_diff[i] / m_k_stiff
            
    # Now calculate the thin filament extensions
    for (i,thin) in enumerate(hs['thin']):
        thin_ext[:,i] = np.diff(thin['bs_x'][::a_bs_per_node])
        
    thin_ext_mean = np.mean(thin_ext, 1)
    thin_ext_sd = np.std(thin_ext, 1)
    
    # And now the thick equivalent
    for (i,thick) in enumerate(hs['thick']):
        thick_ext[:,i] = -np.diff(thick['cb_x'][::m_cbs_per_node])
        
    thick_ext_mean = np.mean(thick_ext, 1)
    thick_ext_sd = np.std(thick_ext, 1)
    
    no_of_rows = 4
    no_of_cols = 1                                      
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             left=0.25,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([7,7])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
            
    ax[0].plot(force_diff)
    ax[0].set_ylabel('Force\nimbalance (nN)')
    ax[0].ticklabel_format(useOffset=False)
    
    ax[1].plot(x_dev)
    ax[1].set_ylabel('Node position\nmismatch (nm)')
    ax[1].set_xlabel('Node number')
    ax[1].ticklabel_format(useOffset=False)

    x = thin['bs_x'][:-2:a_bs_per_node]
    ax[2].plot(x, thin_ext_mean, '-')
    ax[2].fill_between(x, thin_ext_mean-thin_ext_sd,
                       thin_ext_mean+thin_ext_sd,
                       alpha = 0.2)
    ax[2].set_xlim([0, hs['hs_data']['hs_length']])
    ax[2].set_xlabel('Distance along half-sarcomere (nm)')
    ax[2].set_ylabel('Inter-node\nspacing (nm)')
    ax[2].ticklabel_format(useOffset=False)

    x = thick['cb_x'][:-6:m_cbs_per_node]
    ax[3].plot(x, thick_ext_mean, '-')
    ax[3].fill_between(x, thick_ext_mean-thick_ext_sd,
                       thick_ext_mean+thick_ext_sd,
                       alpha = 0.2)
    ax[3].set_xlim([0, hs['hs_data']['hs_length']])
    ax[3].set_xlabel('Distance along half-sarcomere (nm)')
    ax[3].set_ylabel('Inter-node\nspacing (nm)')
    ax[3].ticklabel_format(useOffset=False)
    
    # Work out a file name to save the figure
    file_parts = re.split('/', dump_file_string)
    output_dir = os.path.join(*file_parts[:-3])
    output_file_string = os.path.join(output_dir,
                                      'analysis_%i.png' % analysis_counter)
    # Save the figure
    print('Saving analysis to %s' % output_file_string)
    fig.savefig(output_file_string)
    
            

if __name__ == "__main__":
    
    status_file_strings = [
        '../sim_data/sim_output/1/status_dumps_1_r1/hs_1_time_step_1.json',
        '../sim_data/sim_output/10/status_dumps_1_r1/hs_1_time_step_1.json']    
    
    # Get this directory
    this_dir = os.path.dirname(__file__)
    
    for i, sf in enumerate(status_file_strings):
        
        status_file = os.path.join(this_dir, sf)
        check_force_balance(status_file, i+1)
        