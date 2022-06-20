---
title: Characterization
parent: Structures
has_children: false
nav_order: 6
---

# Characterization

## Overview

Characterization structures simplify
+ comparing models
+ determining standard features including
  + force-pCa curves
  + force-velocity curves

## Details

Characterization structures are written in JSON format and organized into three sections.
+ FiberCpp_exe
+ model
+ characterization

Here is an example from the [thin sweep comparison demo](../../demos/characterization/thin_sweep/thin_sweep.html)

````
{
    "FiberSim_characterization":
    {
        "FiberCpp_exe":
        {
            "relative_to": "this_file",
            "exe_file": "../../../../bin/FiberCpp.exe"
        },
        "model":
        {
            "relative_to": "this_file",
            "model_files": ["models/model_sweep_0.json",
                            "models/model_sweep_10.json",
                            "models/model_sweep_20.json"],
            "options_file": "sim_options.json"
        },
        "characterization":
        [
            {
                "type": "freeform",
                "relative_to": "this_file",
                "protocol_files": ["protocols/twitch_protocol.txt"],
                "sim_folder": "../simulations",
                "m_n": 4,
                "output_image_formats": ["png"],
                "figures_only": "False",
                "trace_figures_on": "True",
                "formatting":
                {
                    "column_titles": [  "Sweep = 0",
                                        "Sweep = 10",
                                        "Sweep = 20"]
                }
            }
        ]
    }
}
````

## FiberCpp_exe

This section points to the binary Cpp file as explained for [batch structures](../batch/batch.html)

## Model

| Key | Comment |
| --- | --- |
| relative_to | root for subsequent paths |
| model_files | an array of [model files](../model/model.html) |
| option file | a base [option file](../options/options.html) |

## Characterization

This section defines the type of characterization to run. Different types of characterization can be run for every model file by providing an array of characterization types.

| Key | Comment |
| --- | --- |
| type | one of, `freeform`, `force_velocity`, `force_pCa` |
| relative_to | root for subsequent paths |
| protocol_files | an array of [protocol files](../protocol/protocol.html), such that each model will be run for every protocol |
| m_n | the number of thick filaments to use in the simulations |
| figures_only | `True` or `False`, setting to true generates the figures from existing simulations and is useful for adjusting the figure options |
| trace_figures_on | `True` or `False`, setting to true generates summary figures for each trial |
