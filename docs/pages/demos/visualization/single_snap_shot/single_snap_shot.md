---
title: Single snap_shot
parent: Visualization
grand_parent: Demos
nav_order: 1
---

# Single snap shot

## Overview

This demo shows how to create a single snap shot using Blender.

## What this demo does

This demo uses the visualization tool from FiberPy to get a single snap shot from a FiberSim simulation.

## Instructions

### Download Blender

First, download Blender at (https://www.blender.org/).

### Getting ready

1. Open an Anaconda Prompt
2. Activate the FiberSim Anaconda Environment by executing:
    ```
    conda activate fibersim
    ```
3. Change directory to `<repo>/code/FiberPy/FiberPy`, where `<repo>` is the directory where you installed FiberSim (e.g. `c:\temp\FiberSim`)

4. Type `python fiberpy.py render_model ../../../demo_files/visualization/single_frame/render_batch.json`

5. A Blender window will open during the process, and close when the rendering is over.

## Viewing the results

+ Use Windows File Explorer to open `<repo>/demo_files/visualization/single_frame/renders`
+ You should see
  + hs_1.png
  
<p align="center">
    <img src="hs_1.png" alt="drawing" width="600"/>
</p>
