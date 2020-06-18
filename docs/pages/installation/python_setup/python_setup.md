---
title: Python Setup
nav_order: 2
parent: Installation
---

# Python Setup
{:.no_toc}

* TOC
{:toc}

FiberSim utilizes the Anaconda platform for managing external python packages. Anaconda is great because it handles all of the external downloads for you, given a special “environment” file. If you would like to learn more about Anaconda and how we in the Campbell Muscle Lab use it, click [this link](http://campbell-muscle-lab.github.io/howtos_Python).

Before running any of the Python files located in the FiberSim repository, we have to give Anaconda an environment file that tells Anaconda which Python packages FiberSim relies on. 

To set up the FiberSim Python environment, follow the instructions located at [this link](https://campbell-muscle-lab.github.io/howtos_Python/pages/anaconda/anaconda.html#using-an-existing-environment), using the `FiberSim\Python\fibersim.yml` file.

## Activating Anaconda Environment

Once you have installed the `FiberSim` environment, to activate your environment:
1. Open an Anaconda Prompt by typing "Anaconda Prompt" into the Windows Start Menu.
2. Type the following in the prompt:
    ```
    conda activate FiberSim
    ```

This gets your Python ready to use all of the packages that FiberPy uses.