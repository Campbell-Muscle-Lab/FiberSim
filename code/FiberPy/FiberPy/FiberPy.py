"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

import sys

from package.modules.display.animate_cb_distributions import animate_cb_distributions as anim_cb
from package.modules.batch.batch import run_batch
from package.modules.batch.batch import run_multiple_batch
from package.modules.output_handler import output_handler as oh
from package.modules.fitting.fit_model import fit_model
from package.modules.visualization import render as viz
from package.modules.visualization.create_movie import create_hs_movie as cm
from package.modules.characterize import characterize_model as characterize
from package.modules.sample import sample_model as sample


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
        fit_model(sys.argv[2])

    if (sys.argv[1] == "render_model"):
        viz.generate_images(sys.argv[2])
    
    if (sys.argv[1] == "movie_with_data"):
        viz.generate_movie_with_data(sys.argv[2])
    
    if (sys.argv[1] == "animate_cb_distribs"):
        if (len(sys.argv) == 5):
            frames = sys.argv[4]
        else:
            frames = []
        anim_cb(
            sys.argv[2], sys.argv[3], frames)

    if (sys.argv[1] == "spatial_visualization"):
        cm(sys.argv[2], sys.argv[3])
    
    if (sys.argv[1] == "characterize"):
        characterize.characterize_model(sys.argv[2])

    if (sys.argv[1] == "run_all_demos"):
        run_multiple_batch(sys.argv[2])
        
    if (sys.argv[1] == "sample"):
        sample.sample_model(sys.argv[2])


if __name__ == "__main__":
    parse_inputs()
