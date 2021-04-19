# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 13:23:35 2021

@author: kscamp3
"""

import numpy as np
import pandas as pd

def create_protocol_file():
    """ Generates a protocol file """

    # Make base_vectors vectors
    n_points = 3000
    time_step = 0.0001

    dt = time_step * np.ones(n_points)
    pCa = 9.0 * np.ones(n_points)
    dhsl = np.zeros(n_points)
    mode = -2 * np.ones(n_points)

    # Implement Ca ste
    pCa[100::] = 5.4

    # Implement force control
    mode[1500:2000] = 10

    # Create pandas dataframe
    df = pd.DataFrame({'dt':dt, 'pCa':pCa, 'dhsl':dhsl, 'mode':mode})
    
    print(df)
    
    # Write to file
    df.to_csv('../sim_input/protocol.txt', sep='\t', index=False)
    
if __name__ == "__main__":
    create_protocol_file()
    
