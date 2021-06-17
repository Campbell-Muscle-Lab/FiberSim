---
title: Job
parent: Structures
has_children: false
nav_order: 10
---

# Job

## Overview

In FiberSim, a `job` defines a single simulation. Jobs are defined in [JSON format](https://en.wikipedia.org/wiki/JSON) and include pointers to additional files.

Here is an example of a job.

````
"job":
[
    {
        "relative_to": "this_file",
        "model_file": "sim_input/model.json",
        "options_file": "sim_input/options.json",
        "protocol_file": "sim_input/pCa45_protocol.txt",
        "results_file": "sim_output/results.txt",
        "output_handler_file": "sim_input/output_handler.json"
    }
]
````

## Details

| Key | Comment |
| ---- | ---- |
| relative_to | root for subsequent paths |
| model_file | the [model_file](../model/model.html) for the simulation |
| options_file | the [options_file](../options/options.html) passed to FiberCpp |
| protocol_file | the [protocol_file](../protocol/protocol.html) passed to FiberCpp |
| results_file | the file to save the [main results](../results/results.html) |
| output_handler | the [output_handler](../output_handler/output_handler.html) which can make figures from the simulation |

