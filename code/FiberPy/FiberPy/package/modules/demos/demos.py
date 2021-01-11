"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

import sys
from package.modules.demos.getting_started.getting_started import parse_inputs as gs_parse_inputs
from package.modules.demos.fitting.fitting import fit_parse_inputs

def parse_inputs(inputs):
    # Get the number of arguments
    no_of_arguments = len(inputs)
    
    if (no_of_arguments > 1):
        if (inputs[2] == "getting_started"):
            print("Getting_started")
            gs_parse_inputs(inputs)
    
        if (inputs[2] == "fitting"):
            print("Fitting")
            fit_parse_inputs(inputs)
                    
 
if __name__ == "__main__":
    parse_inputs()
