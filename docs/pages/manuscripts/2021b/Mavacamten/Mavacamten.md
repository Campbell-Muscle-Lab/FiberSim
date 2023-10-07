---
title: Mavacamten
has_children: False
nav_order: 2
parent: 2021b
grand_parent: Manuscripts
---

# Mavacamten

This page shows how to reproduce a dose-dependent response curve to a myotrope (Mavacamten).

## Getting ready

+ Open an Anaconda Prompt

+ Activate the FiberSim Anaconda Environment by executing:
```
conda activate fibersim
```
+ Change directory to `<FiberSim_dir>/code/FiberPy/FiberPy`, where `<FiberSim_dir>` is the directory where you installed FiberSim.

## Run the simulations

+ Type:
 ```
python FiberPy.py run_batch "../../../manuscripts/2021b/myotrope/batch_mava.json"
 ```

+ You should see text appearing in the terminal window, showing that the simulations are running. When it finishes (this may take ~15 min), you should see something similar to the image below.

![command prompt](command_prompt.PNG)

## Viewing the results

The dose-dependent response is stored in `<FiberSim_dir>/manuscripts/2021b/myotrope/sim_output` 

<img src='mava_dose_response.png' width="50%">

The underlying data are stored in `<FiberSim_dir>/manuscripts/2021b/myotrope/sim_output/analysis.xlsx`

<img src='analysis.PNG' width="40%">

The 10 subfolders from `<FiberSim_dir>/manuscripts/2021b/myotrope/sim_output` contain the simulations results and summary figures for each Mavacamten concentration.

![sim_output](sim_output.PNG)




