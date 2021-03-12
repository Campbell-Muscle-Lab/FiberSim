"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

import sys

from package.modules.batch.batch import run_batch

def parse_inputs():
    # Get the number of arguments
    no_of_arguments = len(sys.argv)

    if (sys.argv[1] == "run_batch"):
        run_batch(sys.argv[2])

    
    # if (no_of_arguments > 0):
    #     if (sys.argv[1] == "demos"):
    #         print("demos")
    #         demos_parse_inputs(sys.argv)

    #     if (sys.argv[1] == "run_my_batch"):
    #         run_my_batch(sys.argv[2])
            
            
                    
 
if __name__ == "__main__":
    parse_inputs()
