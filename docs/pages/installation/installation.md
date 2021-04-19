---
title: Installation
nav_order: 2
---

# Installation

This page discusses how to install FiberSim and the necessary support software.

## Cloning the GitHub repository

The FiberSim code is hosted on GitHub as a repository. You need to download a local copy of this repository on your computer to use FiberSim. To clone the GitHub repository, we recommend using the GitHub Desktop client. You can download the installer for this at [this link](https://desktop.github.com/).

To clone the FiberSim repository, follow the instructions from this video:

[![](installation.PNG)](https://drive.google.com/file/d/1FsMi-QV3wuRo4hMaZE2oW7NQe3Qpn77n/view?usp=sharing)

1. Open GitHub Desktop and click on the "Files" tab.
2. Click on the "Clone repository..." option.
3. Type "campbell-muscle-lab/fibersim" into the search bar.
4. Click on the `Campbell-Muscle-Lab/FiberSim` repository, then choose where you would like the repository to be located on your machine.
5. Click on "Clone".


## Installing Anaconda

Go to <https://www.anaconda.com/products/individual#windows>, and download the “Python 3.8 version.” Accept all default settings. 

## Activating Anaconda Environment 

FiberSim utilizes the Anaconda platform for managing external python packages. Anaconda is great because it handles all of the external downloads for you, given a special “environment” file. If you would like to learn more about Anaconda and how we in the Campbell Muscle Lab use it, click [this link](http://campbell-muscle-lab.github.io/howtos_Python).

Before running any FiberSim simulation, you have to install the FiberSim environment.

1. Open an Anaconda Prompt by typing "Anaconda Prompt" in the Windows Start Menu
2. The command prompt shows the directory you are currently in. Most of the time, the directory is `C:\Users\Name_of_the_user`. You need to change the directory to `<repo>/code/FiberPy/environment`, where `<repo>` is the directory where you cloned the FiberSim repo. To navigate inside the folders, use the "change directory command" (`cd`), as shown here: [![](environment.PNG)](https://drive.google.com/file/d/1L3ANPuRpob6eoh3Vyam4no6U9J6TgpuF/view?usp=sharing)

3. Type the following in the prompt:

```
conda env create -f environment.yml
```

and press `Enter`. Anaconda will handle the download and installation of all dependencies.

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

You are now ready to try the first [demo](../demos/getting_started/getting_started.html)!