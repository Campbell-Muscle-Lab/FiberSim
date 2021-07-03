---
title: Render options
parent: Structures
has_children: false
nav_order: 22
---

# Render options

## Overview

Render options are used for [Visualization](../../demos/visualization/visualization.html).

They tell FiberPy where to find [Blender](https://www.blender.org) and provide some additonal information about the rendering.

## Details

Here is an example of a `render_options file`.

````
{
    "render_options":
    {
        "blender_exe_path": "C:/Program Files/Blender Foundation/Blender 2.92",
        "max_render_distance": 500,
        "max_smooth_distance": 200,
        "background_mode": true
    }
}
````


| Key | Comment |
| ---- | ---- |
| blender_exe_path | the location of the `blender.exe` file |
| max_render_distance | binding sites, myosin heads, and MyBP-C molecules further than from the camera than this distance (in nm) will not be drawn. You can save rendering time by making this number smaller - but your scene will be missing molecules. This might not matter if the molecules would be too far away to be clearly visible anyway. |
| max_smooth_distance | binding sites, myosin heads, and MyBP-C molecules further than from the camera than this distance (in nm) are drawn in less detail. Again, you can save rendering time by making this number smaller. |
| background_mode | true - run without activating the Blender user interface (faster?) or false (activate the interface) |
