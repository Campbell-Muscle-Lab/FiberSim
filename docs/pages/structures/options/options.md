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
	"lambda_jitter" : 20,
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
| max_rate | transition rates calculated faster than this are clipped at `max_rate` - default is 5000 s<sup>-1</sup>|
| x_pos_rel_tol | the tolerance for calculating the x positions in FiberCpp's force-balance calculations - default is 0.001 nm |
| adjacent_bs | the number of binding sites either side of the nearest site to which a myosin head can attach. A value of 0 constrains myosins to attach to the nearest site, while 1 allows a head to attach to any one of 3 sites - default is 0 |
| lambda_jitter | the first myosin crown on each thick filament will be at x = m_lambda + rand * lambda_jitter, where rand is a random number between 0 and 1 - default lambda_jitter is 0|

### Status files

| Key | Comment |
| ---- | ---- |
| relative_to | the root for relative paths |
| status_folder | the folder to store status files |
| time_steps | in the form a:b:c, store status files for time-steps from a to c in increments of b. Thus 1:100:8000 stores snapshots for time-steps 1, 101, 201, 301 to 7901 | 
