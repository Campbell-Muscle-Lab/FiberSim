"""
Single run
@author: Ken Campbell
"""

from pathlib import Path

from .....modules.batch.batch import run_batch
from .....modules.analysis.analysis import analysis as anal

def demo_single_run_with_log():

    # Run a single simulation
    run_batch("package/demo_files/getting_started/single_run_with_log/batch_single_run_with_log.json")
    
    # Results are saved as package/demo_files/temp/single_run/results.txt
   
    # Display the results
    an = anal()
    an.display_data("package/demo_files/temp/single_run/results.txt")

if __name__ == "__main__":
    demo_single_run_with_log()               
    
    
