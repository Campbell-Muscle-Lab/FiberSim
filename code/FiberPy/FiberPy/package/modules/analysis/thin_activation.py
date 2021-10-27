# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 11:38:37 2021

@author: sako231
"""

import os

import numpy as np
import pandas as pd

from package.modules.half_sarcomere import half_sarcomere

def run_thin_analysis(thin_data, batch_file_string):
    """ Creates an excel file containing the RUs states """

    # Pull off the base folder
    base_folder = os.path.dirname(batch_file_string)


    if (thin_data['relative_to'] == 'this_file'):
        data_folder = os.path.join(base_folder,
                                       thin_data['results_folder'])
        output_data_file_string = \
            os.path.join(base_folder, thin_data['output_data_file'])
    else:
        data_folder = thin_data['results_folder']
        output_data_file_string = thin_data['output_data_file']

    
    columns_title = ["Thin Filament ID", "Strand ID (1 or 2)"]
    

    if os.path.isdir(data_folder): # check that the data_folder exists
    
        # Find the results files        
        for file in os.listdir(data_folder):
            
            if file.endswith('.json'): # look for json status files
            
                data_file_string = os.path.join(data_folder, file)
                    
                hs = half_sarcomere.half_sarcomere(data_file_string)
                
                row_counter = 0
                
                for i, thin_fil in enumerate(hs["thin"]): # Loop over all thin filaments
                
                    if i == 0: # Create empty dataframe
                    
                        number_of_ru = int(len(thin_fil["unit_status"])/thin_fil["a_bs_per_node"])
                        
                        for j in range(0, number_of_ru):
                            columns_title.append(f"RU # {j+1}")
                            
                        # Create empty dataframe
    
                        df = pd.DataFrame(columns=columns_title)
                        
                            
                    unit_status = thin_fil["unit_status"][::-1]  # reverse list to start at M-line        
                
                    # First strand contains the odd RUs ! Because we reverted the unit_status list !
                    first_strand = unit_status[1::2]
    
                    # Second strand contains the even RUs
                    second_strand = unit_status[::2]
    
                    # Add thin filament and strand numbers                
                    first_strand.insert(0, i+1) # thin filament number
                    first_strand.insert(1, 1) # strand number
                    
                    second_strand.insert(0, i+1) # thin filament number
                    second_strand.insert(1, 2) # strand number
                   
                    # Now fill two rows of the dataframe
                    
                    
                    df.loc[row_counter,:] = first_strand
                    
                    row_counter += 1
                    
                    df.loc[row_counter,:] = second_strand
                    
                    row_counter += 1
                                    

    # Save the data as an excel file 

    print('Writing data to %s' % output_data_file_string)
    df.to_excel(output_data_file_string,
                engine='openpyxl',
                index=False)