# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 13:23:35 2021

@author: kscamp3
"""

import numpy as np
import pandas as pd

def create_protocol_file():
    """ Generates a protocol file """

    for i in range(0,3):
        # Make base_vectors vectors
        n_points = 12000
        time_step = 0.0001
    
        dt = time_step * np.ones(n_points)
        pCa = 9.0 * np.ones(n_points)
        dhsl = np.zeros(n_points)
        mode = -2 * np.ones(n_points)
    
        # Implement Ca ste
        pCa[500::] = 5.4
    
        # Implement force control
        if (i==0):
            dhsl[0] = 250
        elif (i==1):
            dhsl[2500:7500] = 0.05
        else:
            dhsl[0] = 500
            dhsl[2500:7500] = -0.05
    
        # Create pandas dataframe
        df = pd.DataFrame({'dt':dt, 'pCa':pCa, 'dhsl':dhsl, 'mode':mode})
        
        print(df)
        
        # Write to file
        df.to_csv('../sim_input/protocol_%i.txt' % i, sep='\t', index=False)
    
if __name__ == "__main__":
    create_protocol_file()
    
