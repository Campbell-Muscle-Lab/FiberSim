"""
Single run
@author: Ken Campbell
"""

import json
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path

from .....modules.batch.batch import run_batch
from .....modules.analysis.analysis import analysis as anal

def demo_multiple_run():

    # Run a batch
    batch_file_string = "package/demo_files/getting_started/multiple_run/batch_multiple_run.json"
    run_batch(batch_file_string)
    
    # Make a figure
    f = plt.figure(constrained_layout=True)
    f.set_size_inches([4,4])
    spec2 = gridspec.GridSpec(nrows=1,ncols=1,figure=f)
    ax_0_0 = f.add_subplot(spec2[0,0])
    ax_0_0.set_xlabel('Time (s)')
    ax_0_0.set_ylabel('Force (N m^{-2})')
    
    # Load results back
    with open(batch_file_string) as batch_file:
        json_data = json.load(batch_file)
        batch_data = json_data['FiberSim_batch']
        job_data = batch_data['job']
        for j in job_data:
            print("j_output_folder %s" % j['output_folder'])
            results_file_string = Path(j['output_folder'],'results.txt')
            print("results_file_string: %s" % results_file_string)
            
            fs_data = pd.read_table(results_file_string)
            
            ax_0_0.plot('time', 'force', data=fs_data)
    plt.show()

if __name__ == "__main__":
    demo_multiple_run()               
    
    
