---
title: Visualization
nav_order: 6
parent: Demos
has_children: True
---

# Blender Visualization Pipeline

## Overview

FiberSim can take 'snap-shots' that describe the status of every molecule in the myofibrillar lattice at a given moment. FiberSim's visualization pipeline uses [Blender](https://www.blender.org/) to turn these snap-shots into a 3D scene that can be useful for trouble-shooting or simply watching how the model works.

It's also possible to save a view of the 3D scene as a 2D image. This is called [rendering](https://en.wikipedia.org/wiki/3D_rendering). Rendering successive snap-shots, potentially moving the view-point each time, generates a sequence of images that can ultimately be stitched together to make a movie.

This section provides demos for both of these processes.

## Caveat

Creating a detailed video, like the example below, takes a long time. There are two reasons.

1. You need to write a [status file](../../structures/status_file/status_file.html) to disc for each time-step that you want to show in a movie. This probably takes ~1 s per time-step.

1. You then need to turn each status file into an image. This is very slow. For example, a fast PC will probably take ~10 minutes per thread to render a myofilament lattice with 9 thick filaments.

As an example, if you have 4 threads available on your computer and you want to make a movie with 500 frames (20 seconds at 25 frames per second), this will take ~500 * 10 minutes / 4 = ~1 day to render.

As a result, it's best to plan your videos carefully and run simple tests before launching the final version. It takes a long time to trouble-shoot issues when each attempt requires 1 day.

<iframe width="560" height="315" src="https://www.youtube.com/embed/LMyyscEcL6I" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
