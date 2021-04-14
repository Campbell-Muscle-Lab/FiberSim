---
title: Installation
nav_order: 2
has_children: True
---

# Installation

This page discusses how to install FiberSim and the necessary support software.

## Cloning the GitHub repository

The FiberSim code is hosted on GitHub as a repository. You need to download a local copy of this repository on your computer to use FiberSim. To clone the GitHub repository, we recommend using the GitHub Desktop client. You can download the installer for this at [this link](https://desktop.github.com/).

To clone the FiberSim repo:

1. Open GitHub Desktop and click the dropdown menu for repositories. 
  ![Step 1](step_1.png)
2. Click the "Add" dropdown menu.
  ![Step 2](step_2.png)
3. Click the "Clone repository..." button.
  ![Step 3](step_3.png)
4. Type "FiberSim" into the repository search bar.
  ![Step 4](step_4.png)
5. Click the `campbell-muscle-lab/FiberSim` repository, describe where you would like the repository to be located on your local machine, then click "Clone".
  ![Step 5](step_5.png)

## Installing Anaconda

Go to <https://www.anaconda.com/products/individual#windows>, and download the “Python 3.8 version.” Accept all default settings. 

## Activating Anaconda Environment 

FiberSim utilizes the Anaconda platform for managing external python packages. Anaconda is great because it handles all of the external downloads for you, given a special “environment” file. If you would like to learn more about Anaconda and how we in the Campbell Muscle Lab use it, click [this link](http://campbell-muscle-lab.github.io/howtos_Python).

Before running any FiberSim simulation, you have to install the FiberSim environment. To do so:

1. Start Anaconda Navigator
2. Select the Environments tab (left-hand side)
3. Open a base terminal (left-click on the arrow-head to the right of base(root))
4. Change directory to `<repo>/code/FiberPy/environment/environment.yml`, where `<repo>` is the directory you installed the software (e.g. `c:\temp\FiberSim`)
5. Type the following in the prompt:

```
conda env create -f environment.yml
```

OR

1. Open an Anaconda Prompt by typing "Anaconda Prompt" in the Windows Start Menu
2. Change directory to `<repo>/code/FiberPy/environment/environment.yml`, where `<repo>` is the directory you installed the software (e.g. `c:\temp\FiberSim`)
3. Type the following in the prompt:

```
conda env create -f environment.yml
```
Anaconda will handle the download and installation of all dependencies.

## Using The FiberSim Environment

Each time you want to run FiberSim simulations, you need to launch an Anaconda Prompt to write the command lines. Your first command line should always be to *activate* the FiberSim environment. To do so:

1. Open an Anaconda Prompt
2. Type:

```
conda activate FiberSim
```

You will notice that "base" is now changed to "FiberSim" in the prompt command. 

<p align="center">
<img src="conda_activate.PNG" width="900"/>
</p>

You are now ready to try the first [demo](../Demos/Getting_started/getting_started.html)!