---
title: Simple
parent: Run a Simulation
grand_parent: Demos
nav_order: 1
---

# Running a Simple Simulation
{:.no_toc}

* TOC
{:toc}

## Overview

This demo is a simple demo meant to illustrate the fundamental routine for running FiberSim simulations. We will be simulating a 100 millisecond isometric protocol in pCa 5.5 solution.

<!-- Before you begin, ensure that you have followed the instructions for installing FiberSim [here]() <font color="red">NEED TO UPDATE</font> -->

## Create the Directory Structure for your Simulation

We're going to start off with creating our directory structure for this, and the other demos.

Using Windows File Explorer, navigate to a directory where you would like to store all of your demos, then create a folder called `FiberSim_demos`. Navigate into that folder and  create a directory called `demo_simple`. Your directory structure should now look like the following:

```
C:\temp\FiberSim_demos\
                      |
                      |- demo_simple\
```

We have chosen to setup the `FiberSim_demos` directory in `C:\temp\`, but you can choose to set up this structure anywhere you would like to on your local machine.

Now, copy the FiberSim executable, `FiberSim.exe`, into the `FiberSim_demos` directory that you just made. The structure should now look like this:

```
C:\temp\FiberSim_demos\
                      |
                      |- FiberSim.exe
                      |- demo_simple\
```

## Creating the Instruction Files

For the sake of focusing on running FiberSim and not on the creation of instruction files, we have provided the options, model, protocol, and batch files. Click on the following links and save them under the newly created `demo_simple` directory as `model.json`, `sim_options.json`, and `protocol.txt`.

+ [Model File](model.json)
+ [Options File](sim_options.json)
+ [Protocol File](protocol.txt)

Your directory structure should now look like:

```
C:\temp\FiberSim_demos\
                      |
                      |- FiberSim.exe
                      |- demo_simple\
                      |        |
                      |        |- sim_options.json
                      |        |- model.json
                      |        |- protocol.txt
```

### Creating the Batch File

Now that we have all of the necessary instruction files on the local machine, we can focus on how we actually *run* the FiberSim simulation. FiberSim is ran with 'batch' JSON files that describe the simulation(s) that you are wishing to run. We have included the batch file for running this simulation below, but if you would like to learn more about FiberSim batch files, follow [this link](../../../structures/batch/batch.md).

The structure for our batch file will be the following (with changes to reflect where you have chosen to store your `FiberSim_demos` directory on your machine):

```json
{
    "FiberSim_batch": {
        "job": [
            {
                "model_file_string": "C:\\temp\\FiberSim_demos\\demo_simple\\model.json",
                "options_file_string": "C:\\temp\\FiberSim_demos\\demo_simple\\sim_options.json",
                "protocol_file_string": "C:\\temp\\FiberSim_demos\\demo_simple\\protocol.txt",
                "results_file_string": "C:\\temp\\FiberSim_demos\\demo_simple\\results.txt"
            }
        ]
    }
}
```

Save this in `batch.json` in your local `FiberSim_demos` directory.

## Running FiberSim Batches

To run the FiberSim batch job, do the following from the `FiberSim_demos\demo_simple` directory.

1. Open an Anaconda Prompt.
2. Activate the FiberSim Anaconda Environment by executing:
    ```
    conda activate fibersim
    ```
3. Run the simulation by executing the following in the Anaconda Prompt:
    ```
    python <FIBERSIM_REPO_LOCATION>\Python\modules\FiberSim_utilities.py run_batch batch.json
    ```
    
    Once you hit <kbd>Enter</kbd> on your keyboard, FiberSim will begin setting up and running the simulation you have prepared. FiberSim will print out information pertaining to your simulation to the Anaconda Prompt you have opened. Once it is finished, you can view the results in the newly created `results.txt` file and the dump files located in the `log` folder in the `demo_simple` folder.