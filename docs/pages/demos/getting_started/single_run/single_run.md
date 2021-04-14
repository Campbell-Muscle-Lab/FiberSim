---
title: Single run
has_children: false
parent: Getting started
grand_parent: Demos
nav_order: 1
---

# Single run

## Overview

This demo shows you how to run a simple simulation in FiberSim. 

## What this demo does

This demo simulates a single half-sarcomere that is held isometric and activated in pCa 4.5 solution. Summary results are saved to a data file and plotted as a figure.

## Instructions

### Getting ready

Once you have installed FiberSim:
1. Start Anaconda Navigator
2. Select the Environments tab (left-hand side)
3. Open a FiberSim terminal (left-click on the arrow-head to the right of FiberSim)
4. Change directory to `<repo>/code/FiberPy/FiberPy`, where `<repo>` is the directory you installed the software (e.g. `c:\temp\FiberSim`)

OR

1. Open an Anaconda Prompt
2. Activate the FiberSim Anaconda Environment by executing:
    ```
    conda activate fibersim
    ```
3. Change directory to `<repo>/code/FiberPy/FiberPy`, where `<repo>` is the directory you installed the software (e.g. `c:\temp\FiberSim`)

### Run a simulation

* Type:
 ```
 python FiberPy.py run_batch "../../../demo_files/getting_started/single_run/batch_single_run.json"
 ```

 This command line launches the batch file to run a single job simulation. See [here](../../../FiberPy/structure/structure.html) for more information on the batch structure.


* You should see some text appearing in the terminal window


<p align="center">
<img src="prompt.PNG" width="600"/>
</p>

## Viewing the results

* Use Windows File Explorer to open `<repo>/demo_files/getting_started/single_run/sim_output`
* You should see
  + `results.txt`: contains the raw data from the simulations
  + `summary.png`: summary figure showing the time traces

<br>

<p align="center">
    <img src="summary.png" alt="drawing" width="400"/>
</p>

