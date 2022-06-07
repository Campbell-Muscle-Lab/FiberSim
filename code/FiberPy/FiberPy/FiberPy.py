"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

import sys

from package.modules.batch.batch import run_batch
from package.modules.batch.batch import run_multiple_batch
from package.modules.output_handler import output_handler as oh
from package.modules.fitting import fitting as fit
from package.modules.visualization import render as viz
from package.modules.visualization.create_movie import create_hs_movie as cm
from package.modules.characterize import characterize_model as characterize


def parse_inputs():

    if (sys.argv[1] == "run_batch"):
        if (len(sys.argv)==3):
            run_batch(sys.argv[2])
        else:
            run_batch(sys.argv[2], figures_only=True)

    if (sys.argv[1] == "make_figures"):
        oh.output_handler(sys.argv[2],
                          sim_results_file_string=sys.argv[3])

    if (sys.argv[1] == "fit_model"):
        fit.fitting(sys.argv[2])

    if (sys.argv[1] == "render_model"):
        viz.generate_images(sys.argv[2])
    
    if (sys.argv[1] == "movie_with_data"):
        viz.generate_movie_with_data(sys.argv[2])
    
    # if (sys.argv[1] == "animate_cb_distributions"):
    #     an = anal.analysis()
    #     an.animate_cb_distributions(
    #         sys.argv[2], sys.argv[3])

    if (sys.argv[1] == "spatial_visualization"):
        cm(sys.argv[2], sys.argv[3])
    
    if (sys.argv[1] == "characterize_model"):
        characterize.characterize_model(sys.argv[2])

    if (sys.argv[1] == "run_all_demos"):
        run_multiple_batch(sys.argv[2])


if __name__ == "__main__":
    parse_inputs()
