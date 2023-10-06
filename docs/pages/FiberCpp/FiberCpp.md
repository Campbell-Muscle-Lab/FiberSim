---
layout: default
title: FiberCpp
has_children: true
nav_order: 3
---

# FiberCpp

FiberCpp is the 'core model' of the FiberSim suite. The [code](code/code.html)
+ implements the [calculations](calculations/calculations.html) underlying the simulation
+ is written in C++ for speed
+ is a low-level console application stored as `FiberCpp.exe` in `FiberSim/bin`

FiberCpp is designed to run simulations as quickly and efficiently as possible. It does not have a sophisticated user-interface.

Most people will find it easier to run FiberCpp simulations via [FiberPy](../FiberPy/FiberPy.html). See the [demos](../demos/demos.html) for many examples.

Only those adding new features to the model (or trying to [fix bugs](http://github.com/campbell-muscle-lab/FiberSim/issues)) will need to work directly with FiberCpp. 

## Class diagram

````mermaid
classDiagram
    FiberSim <|-- muscle

    muscle <|-- FiberSim_model
    muscle <|-- FiberSim_options
    muscle <|-- FiberSim_protocol
    muscle <|-- FiberSim_data
    muscle <|-- half_sarcomere
    muscle <|-- series_component

    FiberSim_model <|-- JSON_functions
    FiberSim_model <|-- kinetic_scheme
    
    FiberSim_data <|-- hs_data

    half_sarcomere <|-- thick_filament
    half_sarcomere <|-- thin_filament

    kinetic_scheme <|-- m_state
    
    m_state <|-- transitions
````
