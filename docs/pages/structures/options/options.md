---
title: Options
parent: Structures
has_children: false
nav_order: 18
---

# Options

## Overview

The options file for a [job](../job/job.html) is stored in [JSON format](https://en.wikipedia.org/wiki/JSON) and provides FiberCpp with additional information about how to perform the calculations. It can also be used to tell FiberCpp to save files which describe the status of the model at each time-point in the simulation.

Here is an example.

```` 
{
  "options": {
    "max_rate": 1e4,
    "x_pos_rel_tol": 1e-4,
    "adjacent_bs": 1,
    "status_files": {
      "relative_to": "this_file",
      "status_folder": "../../sim_output/1/hs",
      "time_steps": "1:100:8000"
    }
  }
}
````

## Details

This section explains the components

| Key | Comment |
| ---- | ---- |
| max_rate | transition rates calculated faster than this are clipped at `max_rate` |
| x_pos_rel_tol | the tolerance for calculating the x positions in FiberCpp's force-balance calculations |
| adjacent_bs | the number of binding sites either side of the nearest site to which a myosin head can attach. A value of 0 constrains myosins to attach to the nearest site, while 1 allows a head to attach to any one of 3 sites |

### Status files

| Key | Comment |
| ---- | ---- |
| relative_to | the root for relative paths |
| status_folder | the folder to store status files |
| time_steps | in the form a:b:c, store status files for time-steps from a to c in increments of b. Thus 1:100:8000 stores snapshots for time-steps 1, 101, 201, 301 to 7901 | 
