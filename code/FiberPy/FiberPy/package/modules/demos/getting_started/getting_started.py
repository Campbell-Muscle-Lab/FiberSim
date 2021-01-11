"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

from package.modules.demos.getting_started.single_run.single_run \
    import demo_single_run
from package.modules.demos.getting_started.single_run_with_log.single_run_with_log \
    import demo_single_run_with_log
from package.modules.demos.getting_started.multiple_run.multiple_run \
    import demo_multiple_run

def parse_inputs(inputs):
    # Get the number of arguments
    no_of_arguments = len(inputs)
    
    if (no_of_arguments > 2):
        if (inputs[3] == "single_run"):
            demo_single_run()
        if (inputs[3] == "single_run_with_log"):
            demo_single_run_with_log()            
        if (inputs[3] == "multiple_run"):
            demo_multiple_run()            

if __name__ == "__main__":
    parse_inputs()
