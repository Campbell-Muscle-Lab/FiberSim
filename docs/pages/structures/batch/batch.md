---
title: Batch Files
parent: Structures
nav_order: 2
---

# Batch Files
{:.no_toc}

* TOC
{:toc}

## Overview

This page describes 'batch' files, the main method of running FiberSim simulations.

The batch file is a JSON file that describes either a single simulation or multiple simulations. If you would like to learn more about the JSON file structure, see our writeup [at this link](http://campbell-muscle-lab.github.io/howtos_json).

## Structure Of The Batch File

The structure of the batch file is simple and extensible. There is a `job` JSON array that describes the instruction files and result file of each of the simulations that you would like to run. You can see the basic structure of the batch file below:

```json
{
    "FiberSim_batch": {
        "job": [
            {
                "model_file_string": "FILE\\PATH",
                "options_file_string": "FILE\\PATH",
                "protocol_file_string": "FILE\\PATH",
                "results_file_string": "FILE\\PATH"
            },
            .
            .
            .
            {
                "model_file_string": "FILE\\PATH",
                "options_file_string": "FILE\\PATH",
                "protocol_file_string": "FILE\\PATH",
                "results_file_string": "FILE\\PATH"
            }
        ]
    }
}
```

As you can see:
  + the `job` array contains a single JSON object for each simulation that you wish to run.
  + `job` objects are self-contained, meaning they describe all of the instruction files needed for a single simulation.
  + It's simple to run many jobs at a time. This quite helpful when you would like to simulate families of simulations with the same model file, as is the case when generating tension-pCa curves.

## Running FiberSim Batches

To run a FiberSim batch job, do the following from the directory that contains your batch file:

1. Open an Anaconda Prompt via the Windows Start Menu.
2. Activate the FiberSim Anaconda Environment by executing:
    ```
    conda activate fibersim
    ```
    from within the Anaconda Prompt.
3. Run the simulation by executing the following in the Anaconda Prompt:
    ```
    python <FIBERSIM_REPO_LOCATION>\Python\modules\FiberSim_utilities.py run_batch batch.json
    ```
    
    Once you hit <kbd>Enter</kbd> on your keyboard, FiberSim will begin setting up and running the simulation you have prepared. FiberSim will print out information pertaining to your simulation to the Anaconda Prompt you have opened. Once it is finished, you can view the results in the newly created results file indicated by the `results_file_string` in the batch file you created.