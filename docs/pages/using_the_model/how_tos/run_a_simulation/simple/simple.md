---
title: Simple
parent: Run a Simulation
nav_order: 1
---

# Running a Simple Simulation
{:.no_toc}

* TOC
{:toc}

## Overview

This demo is a simple demo meant to illustrate the fundamental routine for running FiberSim simulations. We will be simulating a 100 millisecond isometric protocol in pCa 5.5 solution.

Before you begin, ensure that you have followed the instructions for installing FiberSim [here]() <font color="red">NEED TO UPDATE</font>

## Create the Directory Structure for your Simulation

We're going to start off with creating our directory structure for this, and the other demos.

Using Windows File Explorer, navigate to a directory where you would like to store all of your demos, then create a folder called `FiberSim_demos`. Navigate into that folder and  create a directory called `demo_a`. Your directory structure should now look like the following:

```
C:\temp\FiberSim_demos\
                      |
                      |- demo_a\
```

We have chosen to setup the `FiberSim_demos` directory in `C:\temp\`, but you can choose to set up this structure anywhere you would like to on your local machine.

Now, copy the FiberSim executable into the `FiberSim_demos` directory that you just made. The structure should now look like this:

```
C:\temp\FiberSim_demos\
                      |
                      |- FiberSim.exe
                      |- demo_a\
```

## Creating the Instruction Files

For the sake of focusing on running FiberSim and not on the creation of instruction files, we have provided the options, model, and protocol files. Click on the following links and save them under the newly created `demo_a` directory as `model.json`, `sim_options.json`, and `protocol.txt`.

+ [Model File](model.json)
+ [Options File](sim_options.json)
+ [Protocol File](protocol_file.txt)

Your directory structure should now look like:

```
C:\temp\FiberSim_demos\
                      |
                      |- FiberSim.exe
                      |- demo_a\
                      |        |
                      |        |- sim_options.json
                      |        |- model.json
                      |        |- protocol.txt
```

### Creating the Batch File

<font color="red">FILL IN HERE</font>

## Running FiberSim Batches

To run the FiberSim batch job, do the following from the `FiberSim_demos` directory.

1. Open an Anaconda Prompt.
2. Activate the FiberSim Anaconda Environment by executing:
    ```
    conda activate fibersim
    ```
3. Run the simulation by executing the following in the Anaconda Prompt:
    ```
    FiberSim.exe demo_a\model.json demo_a\sim_options.json demo_a\protocol.txt demo_a\results.txt
    ```
    The last argument of step 3 tells FiberSim to output the results of the simulation into `demo_a\results.txt`.

    Once you hit <kbd>Enter</kbd> on your keyboard, FiberSim will begin setting up and running the simulation you have prepared. FiberSim will print out information pertaining to your simulation to the Anaconda Prompt you have opened. Once it is finished, you can view the results in the newly created `results.txt` file.