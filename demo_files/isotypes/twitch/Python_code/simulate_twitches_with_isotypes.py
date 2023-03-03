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

    # Variables
    base_setup_file_string = '../base/setup.json'
    base_model_file_string = '../base/base_model.json'

    batch_mode = 'characterize'
    
    generated_model_file_base = '../generated/model'
    generated_setup_file_string = '../generated/generated_setup.json'
    
    FiberSim_code_dir = '../../../../code/fiberpy/fiberpy'

    protocol_file_string = '../protocol/protocol_twitch.txt'

    protocol_module_path = '../../../../code/fiberpy/fiberpy/package/modules/protocols'

    # Different isotypes
    rel_beta = [0, 0.5, 1.0]
    rel_beta_attachment = 0.35

    # Add protocols to the path
    sys.path.append(protocol_module_path)
    import protocols as prot

    # Generate a twitch protocol and write to file    
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
    
    for i in range(len(rel_beta)):

        # Copy the base model
        rep_model = copy.deepcopy(base_model)
        
        # Set the proportion of beta heads
        rep_model['m_parameters']['m_isotype_proportions'] = \
             [1.0-rel_beta[i], rel_beta[i]]

        # Copy the alpha kinetics to beta
        rep_model['m_kinetics'][1] = rep_model['m_kinetics'][0]

        # Change the beta kinetics
        # Attachment
        y = np.asarray(rep_model['m_kinetics'][1]['scheme'][1]['transition'][1]['rate_parameters'],
                  dtype=np.float32)
        y[0] = y[0] * rel_beta_attachment
        rep_model['m_kinetics'][1]['scheme'][1]['transition'][1]['rate_parameters'] = y.tolist()

        # Detachment
        y = np.asarray(rep_model['m_kinetics'][1]['scheme'][2]['transition'][0]['rate_parameters'],
                  dtype=np.float32)
        y[0] = y[0] * rel_beta_attachment
        rep_model['m_kinetics'][1]['scheme'][2]['transition'][0]['rate_parameters'] = y.tolist()

        # Generate a model file name
        rfs = ('%s_%i.json' % (generated_model_file_base, i+1))
        
        # Correct path
        rfs = str(Path(os.path.join(base_dir, rfs)).resolve())
        
        # Check the path exists
        parent_dir = Path(rfs).parent.absolute()
        if not os.path.isdir(parent_dir):
            os.makedirs(parent_dir)

        # Write the model
        with open(rfs, 'w') as f:
            json.dump(rep_model, f, indent=4)
            
        # Convert to absolute path
        rfs = str(Path(os.path.join(base_dir, rfs)).resolve())

        # Add in to array
        rep_model_file_strings.append(rfs)
    
    # Now copy the setup
    generated_setup = copy.deepcopy(base_setup)
    
    # And change the file names to absolute paths
    generated_setup['FiberSim_characterization']['model']['relative_to'] = 'false'
    generated_setup['FiberSim_characterization']['model']['model_files'] = rep_model_file_strings
    generated_setup['FiberSim_characterization']['model']['options_file'] = \
        str(Path(os.path.join(base_dir,
                              generated_setup['FiberSim_characterization']['model']['options_file'])).resolve())
    
    # Loop through the characterizations, changing these paths
    characterize = generated_setup['FiberSim_characterization']['characterization']
    for i,c in enumerate(characterize):
        c['relative_to'] = 'false'
        c['sim_folder'] = str(Path(os.path.join(base_dir,
                                                c['sim_folder'])).resolve())
        p = list()
        for pf in c['protocol_files']:
            p.append(str(Path(os.path.join(base_dir,
                                   pf)).resolve()))
        c['protocol_files'] = p
        generated_setup['FiberSim_characterization']['characterization'][i] = c

    # And finally the setup file
    generated_setup_file_string = str(Path(os.path.join(base_dir,
                                                        generated_setup_file_string)).resolve())
    
    print("\n\n%s\n\n" % generated_setup_file_string)
    
    with open(generated_setup_file_string, 'w') as f:
        json.dump(generated_setup, f, indent=4)
        
    # Generate a command line
    cs = 'pushd \"%s\" & python FiberPy.py %s %s & popd' % \
            (FiberSim_code_dir, batch_mode, generated_setup_file_string)
    
    # And run it
    # os.system(cs)
    
    # Now make figure from the setup structure
    create_figure(generated_setup)
    
def create_figure(setup_structure):
    
    # Variables
    trace_cols = ['r','g','b','m']
   
    # Create figure
    no_of_rows = 6
    fig = plt.figure(constrained_layout=False)
    fig.set_size_inches([3,9])
    spec = gridspec.GridSpec(nrows=no_of_rows,
                              ncols=1,
                              figure=fig,
                              wspace=1,
                              hspace=1)
    ax=[]
    
    # Pull off the sim_data_folder
    sim_folder = setup_structure['FiberSim_characterization']['characterization'][0]['sim_folder']
    top_sim_output_folder = os.path.join(sim_folder, 'sim_output')
    
    output_folders = []
    for item in os.listdir(top_sim_output_folder):
        item = os.path.join(top_sim_output_folder, item)
        if os.path.isdir(item):
            output_folders.append(item)

    # Now loop through pulling up the results files
    for (i,of) in enumerate(output_folders):
        if (i==0):
            for j in range(no_of_rows):
                ax.append(fig.add_subplot(spec[j,0]))

        for item in os.listdir(of):
            if (item == 'rates.json'):
                continue
            else:
                ofs = os.path.join(of, item)
            d = pd.read_csv(ofs, sep='\t')
            
            if (i==0):
                field_names = d.columns
    
            ax[0].plot(d['time'], np.power(10,-d['pCa']))
            ax[1].plot(d['time'], d['hs_length'])
            ax[2].plot(d['time'], d['force'])
            
            for (si,s) in enumerate(['a_pop', 'm_pop', 'c_pop']):
                j = 0
                for fn in field_names:
                    if (fn.startswith(s)):
                        ax[3+si].plot(d['time'], d[fn], '-',
                                   color=trace_cols[j])
                        j=j+1

    # Label
    ax[0].set_ylabel('Ca\n(M)')
    ax[1].set_ylabel('Half-sarcomere\nlength\n(nm)')
    ax[2].set_ylabel('Force\nper\nunit\narea\n(N m^-2)')
    ax[3].set_ylabel('Thin\nfilament')
    ax[4].set_ylabel('Myosin')
    ax[5].set_ylabel('cMYBP-C')
    ax[6].set_xlabel('Time (s)')
    
    
    fig.align_labels()

    # Save the figure
    output_file_string = os.path.join(top_sim_output_folder,
                                      'summary.png')
    fig.savefig(output_file_string, dpi=200, bbox_inches='tight')
    

if __name__ == "__main__":
    demo_twitches_with_isotypes()