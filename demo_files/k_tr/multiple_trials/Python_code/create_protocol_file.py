# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 13:23:35 2021

@author: kscamp3
"""

import numpy as np
import pandas as pd

def create_protocol_file():
    """ Generates a protocol file """

    step_n = 20;
    step_size = 50
    
    pCa_values = np.asarray([5.5, 5.4, 5.25, 4.5])
    
    for i in range(0,(len(pCa_values))):
        # Make base_vectors vectors
        n_points = 10000
        time_step = 0.0001
    
        dt = time_step * np.ones(n_points)
        pCa = 9.0 * np.ones(n_points)
        dhsl = np.zeros(n_points)
        mode = -2 * np.ones(n_points)
    
        # Implement Ca ste
        pCa[100::] = pCa_values[i]
    
        # Implement step
        for j in range(4000, 4000+step_n):
            dhsl[j] = -step_size / step_n
            
        for j in range(4400, 4400 + step_n):
            dhsl[j] = step_size / step_n
        
        for j in range(3900, 4500):
            mode[j] = -1
    
        # Create pandas dataframe
        df = pd.DataFrame({'dt':dt, 'pCa':pCa, 'dhsl':dhsl, 'mode':mode})
        
        print(df)
        
        # Write to file
        df.to_csv('../sim_input/protocol_%i.txt' % i, sep='\t', index=False)
    
if __name__ == "__main__":
    create_protocol_file()
    
