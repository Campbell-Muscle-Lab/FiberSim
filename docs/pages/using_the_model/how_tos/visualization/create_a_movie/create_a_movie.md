---
title: Create a Movie
parent: Visualization
nav_order: 2
---

# Creating an Animation of Your Simulation
{:.no_toc}

* TOC
{:toc}

<font color="red">TO BE UPDATED</font>

Blender outputs single bitmaps (images consisting of pixels) of each time point of the simulation you visualize. This is useful for looking at image stills, but to visualize the dynamics of the model, it is often helpful to make a movie of the simulation. We have included functionality for this in the `animate` submodule of the `util` Python module located in `FiberSim\Python_files\`.

This file is run just as described in ["Running the Visualization Pipeline"](#running-the-visualization-pipeline) except instead of following step 4, change to the `util` directory in `Python_files` and type the following into your command prompt:
```
python animate.py <BITMAP_FILE_ROOT> <OUTPUT_DIRECTORY> <OUTPUT_FILE_NAME> <OUTPUT_ANIMATION_FORMAT> <BITMAP_START_INDEX> <BITMAP_END_INDEX> <PROTOCOL_FILE_PATH>
```
where:
+ `<BITMAP_FILE_ROOT>` is the file path and root for the bitmaps output by Blender. This is the combination of the `output-file-directory` and `output-file-root` optional parameters in the ["Optional Visualization Parameters"](#optional-parameters) section. For example, if `--output-file-directory` was set to `C:\Users\Campbell\scratch` and `--output-file-root` was set to `sample_animation`, `<BITMAP_FILE_ROOT>` would be `C:\Users\Campbell\scratch\sample_animation`.
+ `<OUTPUT_DIRECTORY>` is the directory where the animation will be saved.
+ `<OUTPUT_FILE_NAME>` is the name of the output file.
+ `<OUTPUT_ANIMATION_FORMAT>` is the file format for the output animation of your simulation. See https://imageio.readthedocs.io/en/latest/formats.html#multiple-images for a list of all supported formats (there's a lot). We personally use `.gif`, `.avi`, and `.wmv` a lot.
+ `<BITMAP_START_INDEX>` is the starting index of the sequence of bitmaps you would like to convert into an animation.
+ `<BITMAP_END_INDEX>` is the ending index of the sequence of bitmaps you would like to convert into an animation. Note: This index is not included in the animation.
+ `<PROTOCOL_FILE_PATH>` is the file path to the protocol file you used to make the simulation.

We also have some optional parameters available for the animation utility which are shown here. To add them, simply append `--<PARAMETER_NAME>` to the command line execution:
+ `sampling-freqency` is the frequency with which you would like to show the rendered image. For example, if `sampling-frequency` was set to 5, the animation would show every fifth frame in the sequence of rendered frames between `<BITMAP_START_INDEX>` and `<BITMAP_END_INDEX>`.
+ `frames-per-second` is the number of frames shown per second in the output animation.
