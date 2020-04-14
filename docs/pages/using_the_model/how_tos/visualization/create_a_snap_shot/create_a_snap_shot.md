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

To visualize a model, you must include the `"dump_mode": 1` and the `"log_folder"` simulation options in the options JSON file. These options tell FiberSim to dump the half-sarcomere status files and tells FiberSim where to put those files.

## Creating the Visualization Instructions

To run a visualization, you must pass some instructions to the visualization pipeline in the form of a JSON file.

## Rendering the Snapshot

To render the snapshot, run the following in a command prompt:

```
blender <PATH_TO_FIBERSIM_REPOSITORY>\Python\modules\visualization\generate.py -- -j <VISUALIZATION_INSTRUCTIONS>.json
```