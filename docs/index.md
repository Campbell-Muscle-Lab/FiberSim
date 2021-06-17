---
title: FiberSim
permalink: /
nav_order: 1
---

# FiberSim

FiberSim is software for spatially-explicit modeling of half-sarcomeres. The code tracks the position and status of each myosin head, each binding site on actin, and each molecule of myosin-binding protein C.

<img src="img/FiberSim_render_new.png" width="75%">

You need a Windows PC to run a simulation but you can analyze output from the model on any computer that has an Anaconda-based installation of Python.

## Organization

The main components of the software are:
+ [FiberCpp](pages/FiberCpp/FiberCpp.html) - the core model that implements the calculations underlying the simulations.
  + This software is written in C++ but is currently only compiled for Windows PCs.
+ [FiberPy](pages/FiberPy/FiberPy.html) - accessory code that makes it easier to run different types of simulations, fit models to data, and analyze output.
  + This component is written in Python.

<img src="img/code_structure.png" width="75%">

## Getting started

Check the [demos](pages/demos/demos.html) to see how to:
+ Run simulations
+ Create videos and snap-shots of the model
