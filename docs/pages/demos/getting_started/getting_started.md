---
title: Getting started
has_children: true
parent: Demos
nav_order: 1
---

# Getting started

If this is the first time using FiberSim, make sure you have [installed the software](../../installation/installation.html) before trying the demos. You need a computer running Windows.

## Note

The simulations in the demonstrations have not been fitted to data and are intended to illustrate general mechanisms, as opposed to reproducing specific experimental results.


## Concept


```mermaid

sequenceDiagram
    actor User
    User ->> Hard drive: Creates model file
    Note right of Hard drive: Describes the muscle system
    User ->> Hard drive: Creates protocol file
    Note right of Hard drive: Describes the experiment the user<br>wants to simulate
    User ->> Hard drive: Creates sim_options file
    Note right of Hard drive: Sets advanced options for<br>FiberSim calculations
    User ->> Hard drive: Output handler file (optional)
    Note right of Hard drive: Defines figures to make automatically
    User ->> Hard drive: Creates batch file
    Note right of Hard drive: Groups the<br>model, protocol, options, and output handler files<br>in one place
    User ->> FiberPy: Launches FiberPy with batch file
    FiberPy ->> FiberSim: Runs the simulation<br>defined by the batch file
    Note right of FiberSim: Calculates the behavior<br>of the muscle system
    FiberSim ->> Hard drive: Writes simulation data<br>to a results file
    FiberPy ->>+ Hard drive: If an output handler has been defined,<br>creates figures from results file<br>and save as images
    User ->>+ Hard drive: Views results



```

## Simulations


``` mermaid

sequenceDiagram
    actor User
    User ->> Hard drive: Creates model file
    Note right of Hard drive: Describes the muscle system
    User ->> Hard drive: Creates sim_options file
    Note right of Hard drive: Sets advanced options for<br>FiberSim calculations
    User ->> Hard drive: Creates setup file
    Note right of Hard drive: Defines a set of simulations

    User ->> FiberPy: Launches FiberPy with setup_file
    
    rect rgb(100, 200, 200)
    loop FiberPy creates and organizes<br>all of the files required for<br>set of simulations
        FiberPy ->> Hard drive: protocol files
        FiberPy ->> Hard drive: sim options files
        FiberPy ->> Hard drive: output handlers
    end
    end
    FiberPy ->> Hard drive: batch file
    Note left of Hard drive: Batch file describes a sequence of simulations
    
    rect rgb(200, 100, 100)
    loop Run sequence of simulations
        FiberPy ->> FiberSim: Runs single simulation
        Note right of FiberSim: Calculates the behavior<br>of the muscle system
        FiberSim ->> Hard drive: Writes simulation data<br>to a results file
        FiberPy ->>+ Hard drive: If an output handler has been defined,<br>creates figures from each results file<br>and save as images
    end
    end

    FiberPy ->> Hard drive: Creates figures that summarize the set of simulations
    
    User ->>+ Hard drive: Views results


```