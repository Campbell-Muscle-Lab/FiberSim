---
title: Run your Simulation
parent: Demos
has_children: False
nav_order: 3
---

# Running Your Simulation

## Overview

This section demonstrates how to run your own simulation using the FiberSim model. Before you begin, ensure that you have followed the instructions for installing FiberSim [here](../../installation/installation.md).

## Creating the Directory Structure for your Simulation

We're going to start off with creating a directory structure.

Using Windows File Explorer, navigate to a directory where you would like to store all of your demos, then create a folder called `FiberSim_demos`. Navigate into that folder and create a directory called `my_demo`. Your directory structure should now look like the following:

```
C:\temp\FiberSim_demos\
                      |
                      |- my_demo\
```

We have chosen to setup the `FiberSim_demos` directory in `C:\temp\`, but you can choose to set up this structure anywhere you would like to on your local machine.

Now, copy the FiberSim executable, `FiberSim.exe` (stored in `repo/bin`) , into the `FiberSim_demos` directory that you just made. The structure should now look like this:

```
C:\temp\FiberSim_demos\
                      |
                      |- FiberSim.exe
                      |- my_demo\
```

## Creating the Instruction Files

To run FiberSim from your local machine, you need three instruction files:

+ A model file: this file describe the half-sarcomere model you want to work with (number of myofilaments, kinetic scheme for the crossbridges, ...)
+ An option file: this file describes some computational parameters (such as the required numerical precision) and the log folder path.
+ A protocol file: this file describes the *in silico* protocol you want to simulate, e.g. an isometric twitch, or a ramp stretch.

Three instruction files are provided below as examples. The protocol corresponds to a 100 millisecond isometric contraction in pCa 5.5 solution. Click on the following links and save them under the newly created `my_demo` directory as `model.json`, `sim_options.json`, and `protocol.txt`.

+ [Model File](simple/model.json)
+ [Options File](simple/sim_options.json)
+ [Protocol File](simple/protocol.txt)

The model and option files are JSON files. Those structures are used for FiberSim inputs and outputs, as explained [here](../../structures/structures.md).

Your directory structure should now look like:

```
C:\temp\FiberSim_demos\
                      |
                      |- FiberSim.exe
                      |- my_demo\
                      |        |
                      |        |- sim_options.json
                      |        |- model.json
                      |        |- protocol.txt
```

## Creating the Batch File

Now that we have all of the necessary instruction files on the local machine, we can focus on how we actually *run* the FiberSim simulation. FiberSim is ran with 'batch' JSON files that describe the simulation(s) that you are wishing to run. More information about FiberSim batch files can be found [here](../../../structures/batch/batch.md). We have included the batch file for running the simulation below.

The structure for our batch file will be the following (with changes to reflect where you have chosen to store your `FiberSim_demos` directory on your machine):

```json
{
    "FiberSim_batch": {
        "job": [
            {
                "model_file_string": "C:\\temp\\FiberSim_demos\\my_demo\\model.json",
                "options_file_string": "C:\\temp\\FiberSim_demos\\my_demo\\sim_options.json",
                "protocol_file_string": "C:\\temp\\FiberSim_demos\\my_demo\\protocol.txt",
                "results_file_string": "C:\\temp\\FiberSim_demos\\my_demo\\results"
            }
        ]
    }
}
```

Save this in `batch.json` in your local `FiberSim_demos\my_demo` directory. Your directory structure should now look like:

```
C:\temp\FiberSim_demos\
                      |
                      |- FiberSim.exe
                      |- my_demo\
                      |        |
                      |        |- sim_options.json
                      |        |- model.json
                      |        |- protocol.txt
                      |        |- batch.json
```

## Running FiberSim Batches

To run the FiberSim batch job, do the following from the `FiberSim_demos\my_demo` directory.

1. Open an Anaconda Prompt.
2. Activate the FiberSim Anaconda Environment by executing:
    ```
    conda activate fibersim
    ```
3. Run the simulation by executing the following in the Anaconda Prompt:
    ```
    python <FIBERSIM_REPO_LOCATION>\code\FiberPy\FiberPy\FiberPy.py run_my_batch batch.json
    ```
    
    Once you hit <kbd>Enter</kbd> on your keyboard, FiberSim will begin setting up and running the simulation you have prepared. FiberSim will print out information pertaining to your simulation to the Anaconda Prompt you have opened. Once it is finished, you can view the results in the newly created `results` folder and the dump files located in the `log` folder in the `my_demo` folder.
  