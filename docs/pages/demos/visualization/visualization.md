---
title: Visualization
nav_order: 3
parent: Demos
has_children: True
---

# Blender Visualization Pipeline
{:.no_toc}

* TOC
{:toc}

## Overview

Visualizing mathematical models can be invaluable for debugging the model, gaining insight into the model, or simply creating figures for publications. 

We've developed a pipeline for visualizing the FiberSim model using Blender 2.80 to fit these needs. Blender is an incredibly powerful open source 3D animation suite that supports multiple levels of the 3D model creation process. We created the following animation (inspired by Star Wars!) using the visualization pipeline.

<iframe width="560" height="315" src="https://www.youtube.com/embed/LMyyscEcL6I" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Running the Blender Visualization Pipeline

There are two options for visualizing your FiberSim simulation. You can render a single image as described in the demo [at this link](create_a_snap_shot/create_a_snap_shot.md) or you can render multiple images to create a movie as described in the demo [at this link](create_a_movie/create_a_movie.md).

The general process for visualizing a FiberSim simulation is as follows:

1. Run a FiberSim simulation following guidlines found in the [Running FiberSim](../run_a_simulation/run_a_simulation.md) section.
2. Create a JSON instruction file for your visualization by following the instructions fiound in the [Creating a JSON Instruction File](#creating-a-json-instruction-file) section.
3. Open an Anaconda Prompt and activate the `fibersim` environment.
4. Execute either the `generate.py` file for generating a single snapshot or the python wrapper script, `wrapper.py`, for creating movies of your simulation. Both scripts are found in the `FiberSim\Python\modules\visualization\` directory.

## Creating a JSON Instruction File

The visualization pipeline relies on a user-defined instruction file for the visualization pipeline's parameters. The instruction file utilizes the JSON file format since it is an easily readable and user-friendly format. We recommend downloading and installing Visual Studio Code for instruction file creation. It is entirely free and highlights any errors that may be made in the instruction file.

For more information on the JSON file format and how it is used in FiberSim, see the [JSON Dump Structure Page](../../../structures/json_dump_structure/json_dump_structure.md).

As mentioned previously, there are two approaches to visualization that you can take. You can visualize a single snapshot or you can render multiple snapshots and make them into an animation. The two approaches require instruction files that are similar but have slight differences. For simplicity's sake, we'll start out with describing the structure of the single snapshot instruction file.

### The Single JSON Instruction File

The single JSON instruction file describes the visualization parameters to render a single snapshot of the simulation geometry. The parameters used to describe the visualization are as follows:

+ `dump_file_path`
    + Description: The file path of the half-sarcomere status file that you are wishing to render.
    + Valid Input: Any string describing a FiberSim half-sarcomere status file.
+ `render_file_path`
    + Description: The file path where the rendered image will be saved to.
    + Valid Input: Any string describing a file that ends in `.png`.
+ `render_level`
    + Description: The render level for Blender. `0` is to draw both the thick and thin filaments. `1` is to render only the thick filaments. `2` is to only render the thin filaments.
    + Valid Values: `0`, `1`, or `2`.
+ `render_mode`
    + Description: The render mode for Blender. `0` is to render the model as a biological representation with traditional filament structure while `1` is to render the model as a mechanical representation with springs and nodes. The biological representation is closer to how the sarcomere looks in vivo but mechanical representation is closer to how the model is designed.
    + Valid Values: `0` or `1`.
        + NOTE: The biological render setting is not yet supported.
+ `render_quality`
    + Description: The render quality level you would like to use. You can choose from `medium` or `high`. `medium` uses a very quick render engine that takes less than 5 seconds to render the image at each time point but it doesn't look *quite* as nice as it could. `high` uses an advanced render engine but takes two minutes or more to render the image at each time point.
    + Valid Values: `medium` or `high`
+ `no_render`
    + Description: Turns off the rendering of the model. This is helpful if you would like to just look at one time point of the simulation for debugging purposes.
    + Valid Values: `true` or `false`.
+ `filaments_to_hide`
    + Description: A list of the filaments you would like to "hide" while in Blender and in the rendered screenshots Blender outputs.
    + Valid Values: A list of any number of filament names. The filaments are named as `m_X` and `a_Y` for thick filament with ID #X and thin filament with ID #Y, respectively.
+ `camera`
    + Description: An object that is used to position and orient the camera for rendering images.
    + Example of use:
        ```json
        "camera": {
            "location": {
                "x": 88.119,
                "y": -5.6661,
                "z": 2.1491
            },
            "rotation": {
                "x": 119,
                "y": 5.1,
                "z": 161
            }
        }
        ```
        Where `location` controls the location of the camera using X, Y, and Z coordinates and `rotation` controls the orientation of the camera. The location and orientation can be tricky to guess. See ["Visualization Tips and Tricks"](#visualization-tips-and-tricks) for recommendations on how to place and orient the camera.

#### Sample Single JSON Instruction File

We've included a sample JSON instruction file [here](./sample_single_instruction_file.json).

### The Top-Level JSON Instruction File

The top-level JSON instruction file describes the visualization parameters to render multiple snapshots of the simulation. Much of the top-level JSON structure is the same as described in the [previous section](#the-single-json-instruction-file) with a few minor changes:
  + Instead of `dump_file_path`, we now decribe all of the dump files we would like to visualize with `dump_file_root`, `dump_file_start`, and `dump_file_end`. `dump_file_root` describes the root of the file path to the half-sarcomere dump folder and the half-sarcomere status files you would like to visualize. `dump_file_start` and `dump_file_end` describe the time points that you would like to visualize. For example, if you would like to visualize half-sarcomere status files 1-100 in `C:\temp\log\hs_status`, you would include the following arguments in your top-level JSON file:
      ```json
      "dump_file_root": "C:\\temp\\log\\hs_status\\hs_0_time_step_",
      "dump_file_start": 1,
      "dump_file_end": 101
      ```
  We specify the `dump_file_end` as 101 because the specification is not inclusive (it does not render the final dump file).
  + Instead of `render_file_path`, we now include the `output_file_directory` and `output_file_root` such that the output render will be named `<output_file_root>_<time_step>.png` and will be located in the `output_file_directory`. For example:
        + If we ran our visualization job with `"output_file_directory": "C:\\temp"` and `"output_file_root": "model"`, our renders for the simulation would be `C:\\temp\\model_1.png`, `C:\\temp\\model_2.png`, and so on.
  + We include the file path to the Blender executable in `blender_file_path` like the following:
      ```json
      "blender_file_path": "C:\\Program Files\\Blender Foundation\\Blender\\blender.exe"
      ``` 
  + The `camera` object is now an array of object(s) that include the `location` and `rotation` members. This opens up two options for the array:
      + The array can contain a singular object if you would like for your camera to remain stationary throughout the entirety of the animation. An example of this would be:
          ```json
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
          ```
      + The array can contain as many camera location/rotation objects as there are time points in your visualization. For example, if you are visualizing two time points, you would need two `camera` objects like the following:
          ```json
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
            },
            {
              "location": {
                "x": 1450,
                "y": 2.0,
                "z": 6.0
              },
              "rotation": {
                "x": 90,
                "y": 0,
                "z": 90
              }
            }
          ]
          ```

#### Sample Top-Level JSON Instruction File

We have included a sample top-level JSON instruction file at the [following link](./sample_top_level_instruction_file.json).

## Visualization Tips and Tricks

Creating a cogent visualization is a mixture of science and art, and as such, is typically an iterative process. These are some of the tips and tricks we use to make sure we're getting the most out of every iteration of improving the visualization.

### Preview the Animation.

After simulating the model with FiberSim, it is helpful to visualize a single time step before running a full visualization job. This allows you to find the ideal camera placement and prevents you from wasting time on producing a non-ideal animation.

1. Interact with the 3D model of the half-sarcomere and find a suitable position for the camera (see [this](https://docs.blender.org/manual/en/latest/editors/3dview/navigate/camera_view.html) page for how to position and orient the camera in Blender). 
2. After a suitable position and orientation has been found for the camera, document these in the JSON instruction file created using the instructions found in the [Creating a JSON Instruction File](#creating-a-json-instruction-file) section. 
3. Exit Blender and modify the JSON instruction file to show the simulation at another time step. 
4. Check if the camera is still in an ideal spot and that the animation will be showing you what you want to see. This saves a lot of time by showing you if your view will be blocked later on in the animation. 
5. Once this has been completed, change the JSON instruction file to visualize the entire animation that you would like.

### Hide Filaments That Are Out of the Field of View.
As more filaments are drawn in Blender, the visualization process becomes slower thanks to Blender processing all of the objects that make up the filaments. To mitigate this computational burden, we've added functionality to hide filaments that are out of the field of view of the camera. 
    
To effectively choose which filaments to hide during your animation:
1. Begin by following the first step in this series, [Preview the Animation](#preview-the-animation), and place your camera in an optimal location. 
2. Look through your camera view by pressing <kbd>Numpad0</kbd>. 
3. Continually shift in and out of this view by pressing <kbd>Numpad0</kbd> and select the filaments that are out of the camera's field of view by clicking on one of the objects that make up the filament and hide the filament by clicking the check-box in the upper right hand corner of the screen that corresponds with the filament. This action hides the filament from Blender's main view, the viewport, as well as the render engine, meaning that Blender doesn't have to process all of the objects that make up the filament. 
4. As you carry through with this process, document the filaments that you hid in the `filaments_to_hide` parameter in the JSON instruction file.

## Demos

The following are links to the demos that we have for the visualization pipeline.