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
        "max_threads": 20,
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

If the batch contains more than 1 job, FiberPy can run each job using a separate thread. This means that you can run simulations in parallel which saves wall time.

By default, FiberPy will try to use all but one of the available threads on the machine. Thus, if you have 4 threads available on your PC, FiberPy will use up to three threads at a time (leaving one thread spare to keep your computer "responsive"). If the batch has more than 3 jobs, FiberPy will launch the first 3 together, and then launch subsequent jobs as the initial tasks finish and new threads become available.

If you have a PC with many threads, you might want to keep some threads available for other tasks. You can do this by adding a single key to the batch structure.

| Key | Comment |
| ---- | ---- |
| max_threads | the maximum number of threads FiberPy can use |

Here is an example.

````
    "FiberSim_batch": {
        "max_threads": 20,
        "FiberCpp_exe":
        {
            "relative_to": 0,
            "exe_file": "c:/ken/github/campbellmusclelab/models/fibersim/bin/fibercpp.exe"
        },
       "job":
       [
            {
                "relative_to": "this_file",
                "model_file": "temp\\1\\pCa_650\\model_worker.json",
                "options_file": "sim_input/sim_options.json",
                "protocol_file": "sim_input\\1\\pCa_650\\prot_pCa_650.txt",
                "output_handler_file": "sim_input\\1\\pCa_650\\output_handler_650.json",
                "results_file": "temp\\1\\pCa_650\\results.txt"
            },    
            {
                "relative_to": "this_file",
                "model_file": "temp\\1\\pCa_600\\model_worker.json",
                "options_file": "sim_input/sim_options.json",
                "protocol_file": "sim_input\\1\\pCa_600\\prot_pCa_600.txt",
                "output_handler_file": "sim_input\\1\\pCa_600\\output_handler_600.json",
                "results_file": "temp\\1\\pCa_600\\results.txt"
            }
...
Many jobs
...
            {
                "relative_to": "this_file",
                "model_file": "temp\\1\\pCa_580\\model_worker.json",
                "options_file": "sim_input/sim_options.json",
                "protocol_file": "sim_input\\1\\pCa_580\\prot_pCa_580.txt",
                "output_handler_file": "sim_input\\1\\pCa_580\\output_handler_580.json",
                "results_file": "temp\\1\\pCa_580\\results.txt"
            }
       ]
    }
}
````
