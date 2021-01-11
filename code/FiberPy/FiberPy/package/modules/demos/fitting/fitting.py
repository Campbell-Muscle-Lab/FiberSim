"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

from package.modules.fitting import fitting as fit
from package.modules.batch import batch as bat

def fit_parse_inputs(inputs):
    # Get the number of arguments
    no_of_arguments = len(inputs)
    
    if (no_of_arguments > 2):
        if (inputs[3] == "time_domain_single"):
            demo_time_domain_single()

def demo_time_domain_single():
    """ Demo for fitting a single simulation in the time domain """
    print('Ken was here')
    
    optimization_json_file_string = \
        'package/demo_files/fitting/time_domain_single/optimization.json'
    
    fit_object = fit.fitting(optimization_json_file_string)
    
    
    fit_object.fit_controller()
    
    # print(fit_object.batch_structure)
    
    
    # bat.run_batch(batch_structure = fit_object.batch_structure)
    
    

