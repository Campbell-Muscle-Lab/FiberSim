import os
import sys

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from natsort import natsorted

from pathlib import Path

# Add FiberSim code to path
sys.path.append('../../../../code/fiberpy/fiberpy/package/modules/protocols')


print(sys.path)

import protocols as prot


def create_protocols(release_time_s = 0.17, release_f = []):
    
    this_folder = os.path.abspath(__file__)
    
    for (i, rf) in enumerate(release_f):
        
        pf = prot.create_twitch_protocol(n_points=400)
        
        t = np.cumsum(pf['dt'])
        
        pf.loc[t > release_time_s, 'mode'] = release_f[i]
        
        new_fs = os.path.join(this_folder,
                              '../../protocols/release',
                              'protocol_%i.txt' % (i+1))
        
        prot.write_protocol_to_file(pf, new_fs)
        
    print(pf)
    
if __name__ == "__main__":
    
    create_protocols(release_f = [000, 15000, 30000, 45000, 60000, 75000, 85000, 95000])
        
        
    
    