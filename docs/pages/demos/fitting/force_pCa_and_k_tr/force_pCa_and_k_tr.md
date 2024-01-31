---
layout: default
title: Force-pCa and k_tr
has_children: false
parent: Fitting
grand_parent: Demos
nav_order: 1
---

## Fitting to force-pCa and k_tr

## Overview

This demo shows how to fit a model to force-pCa and k_tr data.

## What this demo does

This demo:

+ uses a [fitting setup](../fitting.html) to iteratively adjust a `base model` seeking the best fit to `target data`.
+ shows how to analyze the fitting process.

## Instructions

If you need help with these step, check the [installation instructions](../../../installation/installation.html).

+ Open an Anaconda prompt
+ Activate the FiberSim environment
+ Change directory to `<FiberSim_repo>/code/FiberPy/FiberPy`
+ Run the command

<b>Note before you run this command that initializing a new fit will over-write the example data stored in your repo. Regenerating it will take a long time - potentially several hours. You may want to read through the rest of this page first.<b>

```text
 python FiberPy.py fit_model "../../../demo_files/fitting/force_pCa_and_k_tr/base/setup.json"
 ```

+ You should see text appearing in the terminal window, showing that the simulations are running. When it finishes (this may take a few minutes), you should see something similar to the image below.

### Viewing the results

All of the results from the simulation are written to files in `<FiberSim_repo>/demo_files/fitting/force_pCa_and_k_tr`

Lots of stuff here.

### How this worked

The model file (found in `<repo>/demo_files/isotypes/k_tr/base/model.json`)defined two isoforms as described [here](../isotypes.html).

The setup file was very similar to that used for the [base k_tr simulation](../../single_trials/k_tr/k_tr.html).

```text
{
  "FiberSim_setup":
  {
    "FiberCpp_exe": {
      "relative_to": "this_file",
      "exe_file": "../../../../bin/FiberCpp.exe"
    },
    "model":
    {
      "relative_to": "this_file",
      "options_file": "sim_options.json",
      "fitting":
      {
        "base_model": "model.json",
        "generated_folder": "../generated",
        "working_folder": "../working",
        "progress_folder": "../progress",
        "Python_objective_call": "../Python_code/return_fit.py",
        "optimizer": "particle_swarm",
        "initial_guess": [0.5, 0.5, 0.5],
        "single_run": "False",
        "adjustments":
        [
            {
                "variable": "m_kinetics",
                "isotype": 1,
                "state": 2,
                "transition": 1,
                "parameter_number": 1,
                "factor_bounds": [-1, 1],
                "factor_mode": "log"
            },
            {
              "variable": "m_kinetics",
              "isotype": 1,
              "state": 3,
              "transition": 1,
              "parameter_number": 1,
              "factor_bounds": [-1, 1],
              "factor_mode": "log"
          },
          {
              "class": "thin_parameters",
              "variable": "a_k_on",
              "output_type": "float",
              "factor_bounds": [0.5, 1.5]
          }
      ]
    }
  },
  "characterization":
  [
        {
            "type": "pCa_length_control",
            "relative_to": "this_file",
            "sim_folder": "../sim_data",
            "m_n": 25,
            "pCa_values": [9, 6.1, 5.8, 5.6, 4.5],
            "sim_duration_s": 1.4,
            "time_step_s": 0.001,
            "pCa_step_up_s": 0.05,
            "k_tr_start_s": 0.8,
            "k_tr_duration_s": 0.02,
            "k_tr_ramp_s": 0.001,
            "k_tr_magnitude_nm": 100,
            "k_tr_fit_time_s": [0.822, 1.39],
            "output_image_formats": [ "png" ],
            "figures_only": "False",
            "trace_figures_on": "False"            
        }
    ]
  }
}
```
