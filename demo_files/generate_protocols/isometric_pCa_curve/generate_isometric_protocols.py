# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:16:13 2022

@author: kscamp3
"""

import os
import sys

import numpy as np

def write_protocols():
    """ Writes a sequence of protocols to a directory and generates a
        batch structure """
        
    # Variables
    protocol_path = 'protocols'
    
    isometric_pCa_values = np.concatenate(([8.0],
                                           np.arange(6.4, 5.0, -0.2),
                                           [4.6]))
    
    # Code
    
    # Add protocols to the path
    sys.path.append('../../../code/fiberpy/fiberpy/package/modules/protocols')
    import protocols
    
    # Create a batch
    FiberSim_batch = dict()
    FiberSim_batch['job'] = []
    
    # Loop through pCa values
    for pCa in isometric_pCa_values:
        
        # Generate the file_string
        prot_file_string = os.path.join(protocol_path,
                                        ('protocol_pCa_%.1f.txt' % (10*pCa)))
        
        # Generate the protocol
        prot = protocols.create_length_control_protocol(time_step = 0.0001,
                                                        n_points = 1000,
                                                        step_pCa_s = 0.001,
                                                        step_pCa = pCa)
        
        # Write protocol to file
        protocols.write_protocol_to_file(prot, prot_file_string)

if __name__ == "__main__":
    write_protocols()