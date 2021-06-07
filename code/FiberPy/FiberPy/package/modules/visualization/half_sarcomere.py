# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:11:28 2020
Modified on Monday Jun 22 by Sarah

@author: kscamp3
"""

import os
import json

import thin_filament as thin
import thick_filament as thick


class half_sarcomere():
    """Half-sarcomere class"""
    
    def __init__(self, json_file_string):
       
        with open(json_file_string) as json_file:
           json_data = json.load(json_file)
 
           # Load hs data
           hs_data = json_data['hs_data']
           self.t_attach_a_node = json_data['titin']['t_attach_a_node']
           self.t_attach_m_node = json_data['titin']['t_attach_m_node']
           self.t_slack_length = json_data['titin']['t_slack_length']
           self.t_k_stiff = json_data['titin']['t_k_stiff']
           self.cb_extensions = hs_data['cb_extensions']
           self.hs_id = hs_data['hs_id']
           self.time = hs_data['time']
           self.hs_length = hs_data['hs_length']
           self.hs_force = hs_data['hs_force']
           self.pCa = hs_data['pCa']
           self.m_nodes_per_thick_filament = \
               hs_data['m_nodes_per_thick_filament']
           self.a_nodes_per_thin_filament = \
               hs_data['a_nodes_per_thin_filament']

           #Load thick filaments
           thick_fil_data = json_data['thick']
           self.thick_fil = []
           for t in thick_fil_data:
               self.thick_fil.append(thick.thick_filament(t))

           #Load thin filaments
           thin_fil_data = json_data['thin']
           self.thin_fil = []
           for t in thin_fil_data:
               self.thin_fil.append(thin.thin_filament(t))
