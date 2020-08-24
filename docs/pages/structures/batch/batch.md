---
title: Batch Files
parent: Structures
nav_order: 1
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
                "results_folder": "FOLDER\\PATH"
            },
            .
            .
            .
            {
                "model_file_string": "FILE\\PATH",
                "options_file_string": "FILE\\PATH",
                "protocol_file_string": "FILE\\PATH",
                "results_folder": "FOLDER\\PATH"
            }
        ]
    }
}
```

As you can see:
  + the `job` array contains a single JSON object for each simulation that you wish to run.
  + `job` objects are self-contained, meaning they describe all of the instruction files and folders needed for a single simulation.
  + It's simple to run many jobs at a time. This quite helpful when you would like to simulate families of simulations with the same model file, as is the case when generating tension-pCa curves.