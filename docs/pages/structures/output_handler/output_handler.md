---
title: Output handler
parent: Structures
has_children: false
nav_order: 19
---

# Output handler

## Overview

A [job](../job/job.html) can contain an optional `output_handler`. This is a file that tells FiberPy to make figures summarizing the simulation when the calculations finish.

Here is an example.

```` 
{
    "templated_images":
    [
        {
            "relative_to": "this_file",
            "template_file_string": "../template/template_summary.json",
            "output_file_string": "../sim_output/summary.png"
        }
    ]
}
````

## Details

This section explains the components

| Key | Comment |
| ---- | ---- |
| relative_to | the root for subsequent file definitions |
| template_file_string | the file defining the figure template (see below) |
| output_file_string | the output file for the figure |

## Templates

FiberPy generates figures from a template defined as a JSON structure.

Here is an example of a template. It will create 3 panels in a single column, with all records showing one field from the [results file](../results/results.html) plotted against time.

Panel 1 shows `pCa`.

Panel 2 shows `hs_length`.

Panel 3 shows:
+ `force`
+ `titin_force`

The templates have many options. For more details, see the [demos](../../demos/demos.html).

````
{
    "x_display":{
        "global_x_field": "time",
        "label": "Time (s)"
    },
    "panels":
    [
        {
            "column": 1,
            "y_info":
            {
                "label": "pCa",
                "ticks": [9, 4.5],
                "series":
                [
                    {
                        "field": "pCa"
                    }
                ]
            }
        },
        {
            "column": 1,
            "y_info":
            {
                "label": "HS length\n(nm)",
                "scaling_type": "close_fit",
                "series":
                [
                    {
                        "field": "hs_length"
                    }
                ]
            }
        },
        {
            "column": 1,
            "y_info":
            {
                "label": "Stress\n(kN m$^{-2}$)",
                "series":
                [
                    {
                        "field": "force",
                        "field_label": "Total",
                        "scaling_factor": 0.001
                    },
                    {
                        "field": "titin_force",
                        "field_label": "Titin",
                        "scaling_factor": 0.001
                    }
                ]
            }
        }
    }
}

| Key | Comment |
| ---- | ---- |
| relative_to | the root for relative paths |
| status_folder | the folder to store status files |
| time_steps | in the form a:b:c, store status files for time-steps from a to c in increments of b. Thus 1:100:8000 stores snapshots for time-steps 1, 101, 201, 301 to 7901 | 
