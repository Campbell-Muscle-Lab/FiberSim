import os
import sys
import json

from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



# Add the validation package

this_path = os.path.dirname(__file__)
code_path = os.path.join(this_path, '../../../../code/FiberPy/fiberpy/package')

code_path = str(Path(code_path).absolute().resolve())

sys.path.append(code_path)

from modules.half_sarcomere import half_sarcomere

def check_force_balance(dump_file_string):
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
    
    print('no_of_thin_filaments: %i' % no_of_thin_filaments)
    print('a_nodes_per_thin_filament: %i' % a_nodes_per_thin_filament)
    print('no_of_thick_filaments: %i' % no_of_thick_filaments)
    print('m_nodes_per_thick_filament: %i' % m_nodes_per_thick_filament)
    print(hs_total_nodes)
    
    force_diff = np.zeros(hs_total_nodes)
    x_dev = np.zeros(hs_total_nodes)
    
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
                      (x_m - x_a - hs['titin']['t_offset'])) + \
                    (hs['titin']['t_sigma'] * \
                     np.exp((x_m - x_a) / hs['titin']['t_L']))
            
            a_row = (thin_id * a_nodes_per_thin_filament) + \
                (hs['titin']['t_attach_a_node'] -1)
            
            force_diff[a_row] = force_diff[a_row] + fd
            force_diff[m_row] = force_diff[m_row] - fd

            if (thin_id==0):
                print('a_row: %i' % a_row)
                
            if (m_f == 0):
                print('m_row: %i' % m_row)
                
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
                
                
                
                
    print('x_m: %6f' % x_m)
    print('x_a: %6f' % x_a)
    

            # Now add in the
            # if (counter == 3403):
            #     print(counter)
            #     print(m_f)
            #     print(i)
            #     print(hs['hs_data']['hs_length'])
            #     print(thick['m_lambda'])
            #     print(x_pos)
            #     print(x)
            #     print(x_neg)
        

    print('max_force_diff: %g\n' % np.amax(np.abs(force_diff)))
    print('max_indices:')
    print(np.argmax(force_diff))
        
    # Calculate deviations
    for i in range(hs_total_nodes):
        if (i < (no_of_thin_filaments * a_nodes_per_thin_filament)):
            x_dev[i] = force_diff[i] / a_k_stiff
        else:
            x_dev[i] = force_diff[i] / m_k_stiff
    
    print('max_x_dev: %g\n' % np.amax(np.abs(x_dev)))        
    
        # node_indices = 
        # for a_n in np.arange(0, len(thin['bs_x']), a_bs_per_node):
        #     print(a_n)
        #     print(thin['bs_x'][a_n])

        #     x = thin['bs_x'][a_n]            
        #     if (a_n == 0):
        #         x_minus = 0
        #         x_plus = thin['bs_x'][a_n + a_bs_per_node]
        #     elif (a_n < (a_nodes_per_thin_filament-a_bs_per_node-1)):
        #         x_minus = thin['bs_x'][a_n - a_bs_per_node]
        #         x_plus = x
        #     else:
        #         x_minus = thin['bs_x'][a_n - a_bs_per_node]
        #         x_plus = thin['bs_x'][a_n + a_bs_per_node]
                
                
        #     fd = (a_k_stiff * (x_plus - x - a_rl)) - \
        #                 (a_k_stiff * (x - x_minus - a_rl))
            
        #     force_diff[counter] = fd
                                 

    no_of_rows = 2
    no_of_cols = 1                                      
    fig = plt.figure(constrained_layout = False)
    spec = gridspec.GridSpec(nrows = no_of_rows,
                             ncols = no_of_cols,
                             figure = fig,
                             wspace = 1,
                             hspace = 1)
    fig.set_size_inches([7,7])
    
    ax = []
    for i in range(no_of_rows):
        for j in range(no_of_cols):
            ax.append(fig.add_subplot(spec[i,j]))
            
    ax[0].plot(force_diff)
    
    ax[1].plot(x_dev)
            
# def return_index(hs, filament, f, n)
#     """ Takes one-based f and n and returns the zero-based
#         index for the node """

#     no_of_thin_filaments = len(hs['thin'])
#     no_of_thick_filaments = len(hs['thick'])
    
#     a_nodes_per_thin_filament = hs['hs_data']['a_nodes_per_thin_filament']
#     m_nodes_per_thick_filament = hs['hs_data']['m_nodes_per_thick_filament']
    
#     a_bs_per_node = hs['thin'][0]['a_bs_per_node']
#     m_cbs_per_node = hs['thick'][0]['m_cbs_per_node']
    
#     thin_nodes = 
            
            
            
            
            
    
    
    
    
    
    
    
    
    

# def run_validation(kinetic_data, batch_file_string):
#     """Entrance function for the FiberSim testing suite"""
    
#     print("Run_validation starting")
    
#     base_folder = os.path.dirname(batch_file_string)
        
#     # Get  model/protocol/option files depending on the validation type
    
#     val_type = kinetic_data["validation_type"]
    
#     if (val_type == 'force_balance'):  # option files is a list for force-balance check
    
#         if isinstance(kinetic_data['options_file'], list): # check if a list of options is provided        
            
#             dump_folder_list = []
#             dump_precision_list = []
        
#             for elmt in kinetic_data['options_file']:
                
#                 if (kinetic_data['relative_to'] == 'this_file'):
                
#                     options_file = os.path.join(base_folder, elmt)
#                     output_folder = os.path.join(base_folder, kinetic_data['output_data_folder'])
                    
#                 else:
                    
#                     options_file = elmt
#                     output_folder = kinetic_data['output_data_folder']                   
                    
#                 with open(options_file, 'r') as f:
#                     opt = json.load(f)
                    
#                 option_dir = os.path.dirname(options_file)
                    
#                 if 'status_files' in opt['options']:
                    
#                     if opt['options']['status_files']['relative_to'] == 'this_file':
#                         dump_folder = os.path.join(option_dir , opt['options']['status_files']['status_folder'])
#                     else:
#                         dump_folder = opt['options']['status_files']['status_folder']
                                    
#                 else:
#                     raise RuntimeError("No dump folder found to run the test")
                    
#                 dump_precision_list.append(opt["options"]["x_pos_rel_tol"])
#                 dump_folder_list.append(dump_folder)
                
#         else:
#             raise RuntimeError("Option file(s) not provided in a list form")          
            
#         # Create output folder if it does not exist
        
#         if not os.path.exists(output_folder):
#             os.makedirs(output_folder) 
                                    
#         print("force_balance check")
        
#         force_balance.compute_force_balance(dump_precision_list, dump_folder_list, output_folder) 
        
#     else:  # it is a kinetic test    
    
#         if (kinetic_data['relative_to'] == 'this_file'):  
            
#             model_file = os.path.join(base_folder, kinetic_data['model_file'])
        
#             protocol_file = os.path.join(base_folder, kinetic_data['protocol_file'])
            
#             options_file = os.path.join(base_folder, kinetic_data['options_file'])
            
#             output_folder = os.path.join(base_folder, kinetic_data['output_data_folder'])
            
#         else:            
#             model_file = kinetic_data['model_file']
#             protocol_file = kinetic_data['protocol_file']
#             options_file = kinetic_data['options_file']
#             output_folder = kinetic_data['output_data_folder']
            
#         # Create output folder if it does not exist

#         if not os.path.exists(output_folder):
#             os.makedirs(output_folder)   
                                       
#         # Get dump_folder from option file
        
#         with open(options_file, 'r') as f:
#             opt = json.load(f)
        
#         option_dir = os.path.dirname(options_file)
                
#         if 'status_files' in opt['options']:
#             if opt['options']['status_files']['relative_to'] == 'this_file':
#                 dump_folder = os.path.join(option_dir, opt['options']['status_files']['status_folder'])
#             else:
#                 dump_folder = opt['options']['status_files']['status_folder']
                
#         else:
#             raise RuntimeError("No dump folder found to run the test")
            
#         # Get number of adjacent bs from option file
        
#         if 'adjacent_bs' in opt['options']:
#             adj_bs = opt["options"]["adjacent_bs"]        
#         else:
#             adj_bs = 0
        
#         # Now run the kinetics test
        
#         val_type = kinetic_data["validation_type"]
        
#         if (val_type == 'm_kinetics'): 
            
#             print("Myosin kinetics check")
            
#             myosin_kinetics.compute_rate(model_file, protocol_file, dump_folder, output_folder, adj_bs)
            
#         if (val_type == 'c_kinetics'): 
                  
#             print("Mybpc kinetics check")
            
#             mybpc_kinetics.compute_rate(model_file, protocol_file, dump_folder, output_folder, adj_bs)
            
#         if (val_type == 'a_kinetics'): 
                    
#             print("Actin kinetics check")
            
#             actin_kinetics.compute_rate(model_file, protocol_file, dump_folder, output_folder)   
            

if __name__ == "__main__":
    
    status_file_string = '../sim_data/sim_output/1/status_dumps_1_1_r1/hs_1_time_step_100.json'    
    
    check_force_balance(status_file_string)
        