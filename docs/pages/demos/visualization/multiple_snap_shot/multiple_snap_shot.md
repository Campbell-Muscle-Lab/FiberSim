---
title: Multiple snap_shot
parent: Visualization
grand_parent: Demos
nav_order: 2
---

# Multiple snap shots

## Overview

This demo shows how to create multiple snap shots using Blender.

## What this demo does

This demo uses the visualization tool from FiberPy to get multiple snap shots from a FiberSim simulation.

## Instructions

### Getting ready

1. Open an Anaconda Prompt
2. Activate the FiberSim Anaconda Environment by executing:
    ```
    conda activate fibersim
    ```
3. Change directory to `<repo>/code/FiberPy/FiberPy`, where `<repo>` is the directory where you installed FiberSim (e.g. `c:\temp\FiberSim`)

4. Type `python fiberpy.py render_model ../../../demo_files/visualization/multi_frame/render_batch.json`

5. Some Blender windows will open during the process, and close when the rendering is over.

## Viewing the results

+ Use Windows File Explorer to open `<repo>/demo_files/visualization/multi_frame/renders`
+ You should see
  + hs_1.png
  + hs_2.png
  + hs_3.png
  + hs_4.png
  + hs_5.png
  
