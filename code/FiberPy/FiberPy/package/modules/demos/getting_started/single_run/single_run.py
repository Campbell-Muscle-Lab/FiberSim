"""
Single run
@author: Ken Campbell
"""

from pathlib import Path

from .....modules.batch.batch import run_batch
from .....modules.analysis.analysis import analysis as anal

def demo_single_run():

    # Run a single simulation
    run_batch("package/demo_files/getting_started/single_run/batch_single_run.json")


if __name__ == "__main__":
    demo_single_run()               
    
    
