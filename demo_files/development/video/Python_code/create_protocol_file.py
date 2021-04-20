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
    step_size = 30
    
    for i in range(0,1):
        # Make base_vectors vectors
        n_points = 25000
        time_step = 0.0001
    
        dt = time_step * np.ones(n_points)
        pCa = 9.0 * np.ones(n_points)
        dhsl = np.zeros(n_points)
        mode = -2 * np.ones(n_points)

    
        # Implement ramp stretch
        for j in range(500,2601):
            dhsl[j] = 0.075
            
        # Implement Ca step
        pCa[3000::] = 5.7

        # Implement step
        for j in range(9500, 9500+step_n+1):
            dhsl[j] = -step_size / step_n
            
        for j in range(9700, 9700 + step_n + 1):
            dhsl[j] = step_size / step_n
        
        for j in range(9400, 9800):
            mode[j] = -1

        # Implement Ca step
        pCa[14500::] = 5.0

        for j in range(19500, 21500):
            mode[j] = 25000

        # Create pandas dataframe
        df = pd.DataFrame({'dt':dt, 'pCa':pCa, 'dhsl':dhsl, 'mode':mode})
        
        print(df)
        
        # Write to file
        df.to_csv('../sim_input/protocol_%i.txt' % i, sep='\t', index=False)
    
if __name__ == "__main__":
    create_protocol_file()
    
