---
title: Create a Snap Shot
parent: Visualization
nav_order: 1
---

# Creating a Snap Shot of Your Simulation
{:.no_toc}

* TOC
{:toc}

## Overview

This demo introduces the workflow for visualizing a single snap shot of a FiberSim simulation with Blender using our homebrewed visualization framework.

This is the first demo in the visualization demo series but it relies on the output of the [Simple Demo](../../run_a_simulation/simple/simple.md). 

The [next demo](../create_a_movie/create_a_movie.md) instructs you on how to create a movie of your simulation.

## Running the Simulation

As stated above, this demo relies on the half-sarcomere dump files of the [Simple Demo](../../run_a_simulation/simple/simple.md). If you have not ran that demo, you must visit that page and follow the instructions stated within to continue.

### Simulation Options Required to Visualize the Model

To visualize a model, you must include the `"dump_mode": 1` and the `"log_folder"` options in the simulation options JSON file. These options tell FiberSim to dump the half-sarcomere status files and tells FiberSim where to put those files.

## The Visualization Instruction File

To run a visualization, you must pass some instructions to the visualization pipeline in the form of a JSON file. The instructions tell the visualization pipeline what geometry to create, how to create the geometry, the location of the camera, and so on. We've included a list of the parameters and short descriptions below. If you would like to read more in depth descriptions about the parameters, see the [Visualization Page](../visualization.md#the-single-json-instruction-file).

  + `dump_file_path` - The file path of the half-sarcomere status file that you are wishing to render.
  + `render_file_path` - The file path where the rendered image will be saved to. This should always end in `.png` to designate that we wish to save in the PNG file format.
  + `render_level` - The render level controls which filaments are drawn in the simulation.
      + Options:
          + `0` - Draw both the thick filaments and thin filaments.
          + `1` - Draw only the thick filaments.
          + `2` - Draw only the thin filaments.
  + `render_mode` - The render mode describes whether you would like to render the animation using the biological representation of the half-sarcomere, using option `0` (not currently supported), or the mechanical representation of the half-sarcomere, using option `1`.
  + `render_quality` - The render quality parameter indicates which render engine you would like to use.
      + Options:
          + `"medium"` - This indicates that you would like to use the very quick but slightly lower quality render engine. Renders with this engine look just fine but are produced in less than a second typically.
          + `"high"` - This indicates that you would like to use the high quality but very slow render engine. Renders with this engine look professional but take several minutes per render. These times stack up significantly when you are producing animations with several thousand frames.
  + `no_render` - A boolean flag that turns off rendering when set to `true`.
  + `filaments_to_hide` - A JSON array that allows for the user to specify certain filaments not to be drawn. The naming convention is `m_XXX` and `a_YYY` to hide thick and thin filament numbers `XXX` and `YYY`, respectively.
  + `draw_mirrored_cb_connections` - The FiberSim model allows for mirrored cross-bridge connections, meaning that a cross-bridge on one end of the lattice is capable of binding to an actin binding site on the other end of the lattice. This is implemented to reduce edge-effects in the simulation. However, this leads to disarray when cross-bridges are drawn attached to an actin binding site on the opposite end of the lattice. To turn this off, set this flag to `false`.
  + `camera` - This is an object describing both the location and orientation of the camera placement for this render. This object has a specific structure that can be seen in the [example JSON file](#example-single-json-instruction-file).
      + `location` - This is the location of the camera in the same units as the dumped half-sarcomere status files. This object has `x`, `y`, and `z` parameters.
          + NOTE: If specifying the location of the camera based off of a manually placed camera in a previously drawn geometry, be sure to take into account the `x_scale` and the `yz_scale` present in the `BlenderParams` class in `generate.py`.
      + `rotation` - This is the rotation of the camera using an Eulerian XYZ rotation scheme. This object has `x`, `y`, and `z` parameters in degrees. 

### Example Single JSON Instruction File

For simplicity's sake, we have provided an example render file for this demo. The JSON instruction file would look similar to this:

```json
{
  "dump_file_path": "C:\\temp\\log\\hs_status\\hs_0_time_step_1.json",
  "render_file_path": "C:\\temp\\output\\renders\\model_1.png",
  "render_level": 0,
  "render_mode": 1,
  "render_quality": "medium",
  "no_render": false,
  "filaments_to_hide": [],
  "draw_mirrored_cb_connections": false,
  "camera": {
    "location": {
      "x": 1573.9025,
      "y": 0.9409090909090908,
      "z": 5.797357954545454
    },
    "rotation": {
      "x": 79.2145,
      "y": 0.0,
      "z": 90.0
    }
  }
}
```

## Rendering the Snapshot

To render the snapshot, run the following in a command prompt:

```
<PATH_TO_BLENDER_EXE> <PATH_TO_FIBERSIM_REPOSITORY>\Python\modules\visualization\generate.py -- -s <VISUALIZATION_INSTRUCTIONS>.json
```

The `-s` flag indicates that this is a single JSON instruction file to produce a singular render.