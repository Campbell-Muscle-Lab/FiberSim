---
title: Visualization
nav_order: 3
parent: How To's
grand_parent: Using the Model
has_children: True
---

# Blender Visualization Pipeline
{:.no_toc}

* TOC
{:toc}

## Overview

Visualizing mathematical models can be invaluable for debugging the model, gaining insight into the model, or simply creating figures for publications. 

We've developed a pipeline for visualizing the FiberSim model using Blender 2.80 to fit these needs. Blender is an incredibly powerful open source 3D animation suite that supports multiple levels of the 3D model creation process.

<font color="red">To see how to install Blender, visit the ["Installing Blender"](../installation/installing_blender/installing_blender.md) section.</font>

## Running the Blender Visualization Pipeline

Visualizing the model is quite simple. To do so, follow these steps:

1. Run a FiberSim simulation following guidlines found in the <font color="red">[Running FiberSim](../running_fibersim/running_fibersim.md)</font> section.
2. Create a JSON instruction file for your visualization by following the instructions found in the [Creating a JSON Instruction File](#creating-a-json-instruction-file) section.
3. Open either command line or terminal.
    + For Windows users, you'll type `cmd` into the start menu.
    + For Mac and Linux users, type `terminal` into your operating system search bar.
4. Change your directory to the `visualization` directory in the `Python\modules` directory using the `cd` command.
    + Ex: `cd C:\Users\Campbell\Dylan\GitHub\Models\FiberSim\Python\modules\visualization`
5. Execute the python script that tells Blender how to draw the model. The script is, `generate.py` and is ran by typing the following into the command prompt:
    ```
    blender --python generate.py -- -j <JSON_INSTRUCTION_FILE>
    ```
    where `<JSON_INSTRUCTION_FILE>` is the file path to the instruction file you created in step 2.

## Creating a JSON Instruction File
The visualization pipeline relies on a user-defined instruction file for the visualization pipeline's parameters. The instruction file utilizes the JSON file format since it is an easily readable and user-friendly format. We recommend downloading and installing Visual Studio Code for instruction file creation. It is entirely free and highlights any errors that may be made in the instruction file.

To create an instruction file:
1. Open your text editor of choice. Again, our editor of choice for JSON files is Visual Studio Code.
2. Open a new file and name it whatever you would like, including the `.json` file extension.
3. Include all of the required parameters and any of the optional parameters you would like to. All of which are found below.

### Required Parameters
These are the parameters that *must* be in the JSON instruction file that you create.
+ `dump_file_root`
    + Description: The file root of the simulation dump files. This includes the file path of the dump files and the name of the file up until the number of the time step of the dump file. For example, if your simulation dumped files in `C:\temp\log\hs_status\` and the dump files were `hs_0_time_step_XXX.json` where `XXX` was the time step number, then `dump_file_root` would be `"C:\\temp\\dump\\hs_0_time_step_"`.
    + Valid Input: Any string pointing to a directory where FiberSim has output the dump files of a simulation.
+ `dump_file_start`
    + Description: The number indicating where the time point with which you would like the model visualization to start.
    + Valid Input: Any integer that corresponds to a dump file at that time point.
+ `dump_file_end`
    + Description: The number indicating the time point with which you would like the model visualization to stop. Note: this is non-inclusive meaning that Blender will not render the time point indicated here.
+ <font color="red">`camera` - need to update </font>
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

It is important to note that the above are the only *required* parameters. However, the visualization pipeline we have created supports multiple options on how to create the visualization. See the ["Optional Parameters"](#optional-parameters) section to see a list of the name and function of these.

### Optional Parameters
These parameters are optional and thus have default values if none are specified.

+ `output_file_directory`
    + Description: Where the output files of the render will be saved.
    + Default Value: `FiberSim\Python_files\visualize\blender_renders\`
    + Valid Values: Any string describing a directory path. Note: If the backslash is used to specify a path, you must use two backslashes for each backslash in the path.
+ `output_file_root`
    + Description: The base name of all of the rendered images that are output.
        + Ex: if this value is `img`, all rendered images would be saved as img_1, img_2, etc.
    + Default Value: `model`
    + Valid Values: Any string without forward- or back-slashes.
+ `render_level`
    + Description: The render level for Blender. '0' is to render the entire animation. '1' is to render only the thick filaments. '2' is to only render the thin filaments.
    + Default Value: 0
    + Valid Values: 0, 1, or 2.
+ `render_mode`
    + Description: The render mode for Blender. '0' is to render the model as a biological representation with traditional filament structure while '1' is to render the model as a mechanical representation with springs and nodes. The biological representation is closer to how the sarcomere looks in vivo but mechanical representation is closer to how the model is designed.
    + Default Value: 1
    + Valid Values: 0 or 1.
+ `render_quality`
    + Description: The render quality level you would like to use. You can choose from `medium` or `high`. `medium` uses a very quick render engine that takes less than 5 seconds to render the image at each time point but it doesn't look *quite* as nice as it could. `high` uses an advanced render engine but takes two minutes or more to render the image at each time point.
    + Default Value: `medium`
    + Valid Values: `medium` or `high`
+ `no_render`
    + Description: Turns off the rendering of the model. This is helpful if you would like to just look at one time point of the simulation for debugging purposes.
    + Default Value: `false`
    + Valid Values: `true` or `false`.
+ `filaments_to_hide`
    + Description: A list of the filaments you would like to "hide" while in Blender and in the rendered screenshots Blender outputs.
    + Default Value: [ ] (an empty list)
    + Valid Values: A list of any number of filament names. The filaments are named as `m_X` and `a_X` for thick filament with ID #X and thin filament with ID #X, respectively.

Doing the above steps outputs bitmaps into the folder specified by `output_file_directory`. To make these images into a movie, we have included instructions to use the included utility in the ["Creating a Movie of Your Simulation"](#creating-a-movie-of-your-simulation) section.

### Sample JSON Instruction File

We've included a sample JSON instruction file [here](./sample_instruction_file.json).

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