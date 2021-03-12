"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

import sys

from package.modules.batch.batch import run_batch
from package.modules.output_handler import output_handler as oh

def parse_inputs():

    if (sys.argv[1] == "run_batch"):
        run_batch(sys.argv[2])

    if (sys.argv[1] == "make_figures"):
        oh.output_handler(sys.argv[2],
                          sim_results_file_string=sys.argv[3])

if __name__ == "__main__":
    parse_inputs()
