---
title: FiberCpp
has_children: true
nav_order: 3
---

# FiberCpp

FiberCpp is the 'core model' of the FiberSim suite. The [code](code/code.html)
+ implements the [calculations](calculations/calculations.html) underlying the simulation
+ is written in C++ for speed
+ is a low-level console application stored as `FiberSim.exe` in `repo/bin`

FiberCpp is designed to run simulations as quickly and efficiently as possible. It does not have a sophisticated user-interface.

Most people will find it easier to run FiberCpp simulations via [FiberPy](../FiberPy/FiberPy.html). See the [demos](../demos/demos.html) for many examples.

Only those adding new features to the model (or trying to [fix bugs](http://github.com/campbell-muscle-lab/FiberSim/issues)) will need to work directly with FiberCpp. 