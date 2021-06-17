---
title: Batch
parent: Structures
has_children: false
nav_order: 3
---

# Batch

## Overview

FiberPy launches simulations using a batch structure stored in [JSON format](https://en.wikipedia.org/wiki/JSON).

The batch structure contains information about the FiberCpp core model and at least 1 [job](../job/job.html). Each job is a separate simulation.

Here is a simple example from the [isometric_activation_demo](../../demos/getting_started/isometric_activation/isometric_activation.html).

````
{
    "FiberSim_batch": {
        "FiberCpp_exe":
        {
            "relative_to": "this_file",
            "exe_file": "../../../bin/FiberCpp.exe"
        },
        "job":[
            {
                "relative_to": "this_file",
                "model_file": "sim_input/model.json",
                "options_file": "sim_input/options.json",
                "protocol_file": "sim_input/pCa45_protocol.txt",
                "results_file": "sim_output/results.txt",
                "output_handler_file": "sim_input/output_handler.json"
            }
        ]
    }
}
````

## Details

This section explains the components

### FiberCpp_exe

This section points to the binary (*.exe) of the FiberCpp core model.

| Key | Comment |
| ---- | ---- |
| relative_to | root for subsequent paths |
| exe_file | the path to the FiberCpp.exe file |

### Jobs

The [job page](../job/job.html) provides information about the job structure.

## Multiple jobs

The batch structure can contain more than 1 job. Here is an example with 2.

````
{
    "FiberSim_batch": {
        "FiberCpp_exe":
        {
            "relative_to": "this_file",
            "exe_file": "../../../bin/FiberCpp.exe"
        },
        "job":[
            {
                "relative_to": "this_file",
                "model_file": "sim_input/model.json",
                "options_file": "sim_input/options.json",
                "protocol_file": "sim_input/pCa90_protocol.txt",
                "results_file": "sim_output/results_1.txt",
                "output_handler_file": "sim_input/output_handler_1.json"
            },
            {
                "relative_to": "this_file",
                "model_file": "sim_input/model.json",
                "options_file": "sim_input/options.json",
                "protocol_file": "sim_input/pCa45_protocol.txt",
                "results_file": "sim_output/results_2.txt",
                "output_handler_file": "sim_input/output_handler_2.json"
            }
        ]
    }
}
````

The [model](../model/model.html) and the [options](../options/options.html) are identical for both jobs but they use different [protocols](../protocol/protocol.html), save the [results](../results/results.html) to different files, and use different [output_handlers](../output_handler/output_handler.html)

In this case, the protocol filenames show that the two jobs define simulations for pCa 9.0 and pCa 4.5 respectively. If additional jobs defining simulations for intermediate pCa values were added, the batch structure could be used to generate a tension-pCa curve. See the [pCa curve demos](../../demos/pCa_curves/pCa_curves.html) for examples.

### Parallel processing

If the batch contains more than 1 job, FiberPy will attempt to run each job in parallel as a separate thread using all but one of the available processors on the machine.

Thus, if the PC has 4 processors, you can run 3 jobs in almost the same time as you can run one.

If the batch has more than 3 jobs, FiberPy will run the jobs in sequence based on the availability of processors, in this case starting the first 3 jobs simultaneously, and then launching the 4th job once the first simulation has finished.

Note that your computer will other programs more slowly if it is using most of its processors to run FiberSim calculations.

