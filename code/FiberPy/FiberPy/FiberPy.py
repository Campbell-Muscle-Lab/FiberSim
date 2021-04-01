"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

import sys

from package.modules.batch.batch import run_batch
from package.modules.output_handler import output_handler as oh
from package.modules.fitting import fitting as fit
from package.modules.visualization import render as viz

def parse_inputs():

    if (sys.argv[1] == "run_batch"):
        run_batch(sys.argv[2])

    if (sys.argv[1] == "make_figures"):
        oh.output_handler(sys.argv[2],
                          sim_results_file_string=sys.argv[3])

    if (sys.argv[1] == "fit_model"):
        fit.fitting(sys.argv[2])

    if (sys.argv[1] == "render_model"):
        viz.generate_images(sys.argv[2])


if __name__ == "__main__":
    parse_inputs()
