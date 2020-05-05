---
title: Create a Movie
parent: Visualization
nav_order: 2
---

# Creating an Animation of Your Simulation
{:.no_toc}

* TOC
{:toc}


## The Visualization Wrapper

In our [previous demo](../create_a_snap_shot/create_a_snap_shot.md), we create a single bitmap of our simulation. This is useful for looking at a single time point in the simulation, but to visualize the dynamics of the model, it is often helpful to make a movie of the simulation. We have included functionality for this in the `wrapper.py` file found in `FiberSim/Python/modules/visualization/`.

## The Top-Level Instruction File

To run this file, you must generate a top-level JSON file that will describe the entire visualization job. This top-level JSON will then be split into multiple JSON files that will translate to one rendered image per instruction file. This allows minor tweaks to be made in camera position or rotation if so desired.

### Top-Level Structure

Much of the top-level JSON structure is the same as described in the [previous demo](../create_a_snap_shot/create_a_snap_shot.md#the-visualization-instruction-file) with a few minor changes that are summarized [here](../visualization.md#the-top-level-json-instruction-file).

### Our Top-Level Instruction File

We have included the top-level JSON file below:

```json
{
  "dump_file_root": "C:\\temp\\log\\hs_status\\hs_0_time_step_",
  "dump_file_start": 1,
  "dump_file_end": 101,
  "blender_file_path": "C:\\Program Files\\Blender Foundation\\Blender\\blender.exe",
  "output_file_directory": "output\\renders",
  "output_file_root": "model",
  "render_level": 0,
  "render_mode": 1,
  "render_quality": "medium",
  "no_render": false,
  "filaments_to_hide": [],
  "camera": [
    {
      "location": {
        "x": 1500,
        "y": 1.5,
        "z": 5.8
      },
      "rotation": {
        "x": 90,
        "y": 0,
        "z": 90
      }
    }
  ]
}
```

## Rendering the Movie

To render the movie:

1. Open an Anaconda Prompt.
2. Activate the `fibersim` environment.
3. Navigate to the directory where the top-level JSON file is located.
4. Execute the following:
    ```
    python <PATH_TO_FIBERSIM_REPOSITORY>\Python\modules\visualization\wrapper.py -- -m <TOP_LEVEL_JSON>.json
    ```

    Where the `-m` flag indicates that the job will be outputting multiple renders.