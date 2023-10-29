---
layout: default
title: Variable half-sarcomeres
has_children: false
parent: Myofibrils
grand_parent: Demos
nav_order: 3
---

# Half-sarcomeres with variable properties

## Overview

This demo shows how to simulate multiple half-sarcomeres that have variable properties.

## What this demo does

This demo:

+ Builds on the [multiple half-sarcomeres with series compliance](../multiple_hs_with_sec/multiple_hs_with_sec.html) demo but reduces the attachment rate for myosin in a half-sarcomere in the middle of the myofibril
+ Plots summaries of the simulation

## Instructions

If you need help with these step, check the [installation instructions](../../../installation/installation.html).

+ Open an Anaconda prompt
+ Activate the FiberSim environment
+ Change directory to `<FiberSim_repo>/code/FiberPy/FiberPy`
+ Run the command
```
 python FiberPy.py characterize "../../../demo_files/myofibrils/variable_hs/base/setup.json"
```

+ You should see text appearing in the terminal window, showing that the simulations are running. It may take a few minutes to finish.

### Viewing the results

All of the results from the simulation are written to files in `<FiberSim_repo>/demo_files/myofibrils/variable_hs/sim_data/sim_output`

The file `superposed_traces.png` shows pCa, length, force per cross-sectional area (stress), and thick and thin filament properties for the first half-sarcomere in the myofibril plotted against time. Note the complex time-course of relaxation.

<img src="images/superposed_traces.png" width="50%">

FiberPy also generated a customized figure for this simulation showing the force and length of each of the 5 half-sarcomeres in the myofibril. This file is named `summary.png` and is saved to the normal `sim_output` folder.

<img src="images/summary.png" width="50%">

### How this worked

The demo was identical to that shown in [multiple half-sarcomeres with series compliance](../multiple_hs_with_sec/multiple_hs_with_sec.html) except that the number of half-sarcomeres was reduced to 5 and the rate of the attachment power stroke was reduced in the 3rd half-sarcomere in the myofibril.

The top of the model file (found at `<repo>/demo_files/myofibrils/variable_hs/base/model.json`) was:

```
"muscle": {
    "no_of_half_sarcomeres": 5,
    "no_of_myofibrils": 1,
    "sc_k_stiff": 3000,
    "initial_hs_length": 1100,
    "prop_fibrosis": 0.0,
    "prop_myofilaments": 0.5,
    "m_filament_density": 0.407e15
  },
  "half_sarcomere_variation":
  [
    {
      "variable": "m_kinetics_isotype_1_state_2_transition_2_parameter_1",
      "multiplier": [1, 1, 0.6, 1, 1]
    }
  ]
```

Note that `no_of_half_sarcomeres` was set to 5.

A new section titled `half_sarcomere_variation` defined an array of manipulations. In this example, the variable name points to:
+ myosin isotype 1
+ state 2
+ transition 2
+ parameter 1

The m_kinetic scheme is defined further down the model file. The keys above correspond to a parameter value of 200.

```
"m_kinetics": [
  {
    "no_of_states": 4,
    "max_no_of_transitions": 2,
    "scheme": [
      {
        "number": 1,
        "type": "S",
        "extension": 0,
        "transition": [
          {
            "new_state": 2,
            "rate_type": "force_and_mybpc_dependent",
            "rate_parameters": [ 20, 200, 1, 1, 0.01, 0.01, 1, 1]
          }
        ]
      },
      {
        "number": 2,
        "type": "D",
        "extension": 0,
        "transition": [
          {
            "new_state": 1,
            "rate_type": "constant",
            "rate_parameters": [ 100 ]
          },
          {
            "new_state": 3,
            "rate_type": "gaussian_hsl",
            "rate_parameters": [ 200 ]
          }
        ]
      },
      {
        "number": 3,
        "type": "A",
        "extension": 0.0,
        "transition": [
          {
            "new_state": 2,
            "rate_type": "poly",
            "rate_parameters": [ 500, 10, 4]
          },
          {
              "new_state": 4,
              "rate_type": "constant",
              "rate_parameters": [ 75 ]
            }
        ]
      },
      {
          "number": 4,
          "type": "A",
          "extension": 5.0,
          "transition": [
            {
              "new_state": 2,
              "rate_type": "exp_wall",
              "rate_parameters": [ 75, 1, 8, 7]
            }
          ]
        }
    ]
```
The `characterization` was the same as for the [single half-sarcomere simulation](../hs_with_sec/hs_with_sec.html) except that a `post_sim_Python_call` was added. This points to a a standard Python file that is called after the simulation has finished. In this case, the Python made the `summary.png` figure showing the responses of the different half-sarcomeres.

```
{
  "FiberSim_setup":
  {
    "FiberCpp_exe": {
      "relative_to": "this_file",
      "exe_file": "../../../../bin/FiberCpp.exe"
    },
    "model": {
      "relative_to": "this_file",
      "options_file": "sim_options.json",
      "model_files": ["model.json"]
    },
     "characterization": [
        {
            "type": "pCa_length_control",
            "relative_to": "this_file",
            "sim_folder": "../sim_data",
            "m_n": 9,
            "pCa_values": [4.5],
            "sim_duration_s": 1.0,
            "time_step_s": 0.001,
            "pCa_step_up_s": 0.1,
            "pCa_step_down_s": 0.7,
            "output_image_formats": [ "png" ],
            "figures_only": "False",
            "trace_figures_on": "False",
            "post_sim_Python_call": "../Python_code/multiple_hs_summary.py"
        }
    ]
  }
}
```

FiberCpp sets the value of this parameter in each half-sarcomere to the base value (200) multiplied by the corresponding value in the `multiplier` array. Thus:
+ half-sarcomere 1, 1 * 200
+ half-sarcomere 2, 1 * 200 = 200
+ half-sarcomere 3, 0.6 * 200 = 120
+ half-sarcomere 4, 1 * 200 = 200
+ half-sarcomere 5, 1 * 200 = 200
