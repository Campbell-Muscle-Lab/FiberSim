---
title: FiberCpp
has_children: true
nav_order: 3
---

# FiberCpp

FiberCpp is the 'core model' of the FiberSim suite. The code
+ implements the calculations underlying the simulation
+ is written in C++ for speed
+ is a low-level console application stored as `FiberSim.exe` in `repo/bin`

FiberCpp is designed to run simulations as quickly and efficiently as possible. It does not have a sophisticated user-interface.

Most people will find it easier to run FiberCpp simulations via [FiberPy](../FiberPy/FiberPy.html). See the [demos](../demos/demos.html) for many examples.

Note that users can implement different types of kinetic scheme for both myosin and myosin binding protein-C using the [kinetic_scheme](../structures/model/kinetic_scheme/kinetic_scheme.html) in the [model_file](../structures/model/model.html). Only those adding new features to the model (or trying to fix bugs!) will need to work directly with the FiberCpp source code. 