---
layout: default
title: Fitting
has_children: true
parent: Demos
nav_order: 8
---

## Fitting

In FiberSim, _fitting_ a model means trying to find the set of parameter values that yields the optimal match to some _target data_.

In principle, fitting sounds relatively easy - few readers of this documentation would struggle to find a straight line to data, which implies finding the best combination of `m` and `c` for `y = m x + c`. In practice, fitting FiberSim models is much more challenging, and sometimes feels more like an art than a science.

There are several reasons:

+ FiberSim models take a long time (typically tens of seconds) to run, so testing each potential combination of parameters is time-consuming.
+ There are _lots_ of parameters one could adjust - how does one decide which parameters to try?
+ FiberSim is stochastic, so each simulation has some (perhaps random?) uncertainty associated with it.
  + The uncertainty can be reduced by averaging over more filaments, but this adds to the computational time.
+ FiberSim is complicated so each test involves generating, or at least organizing, different combinations of model, protocol, and options files.
+ How should the fit between the simulation and the target data be evaluated?

This set of demonstrations provides some examples. If you want to fit your own models, it's likely that you are involved some sort of active research. Please reach out if you need more help with a specific idea.

### Overview

The basic sequence is as follows.

<!---
sequenceDiagram
    actor User 
    User->>+FiberPy: Provides setup, protocol, base model files<br/>to run simulations
    User->>+FiberPy: Provides target data
    User->>+FiberPy: Provides Python code to quantify fit between<br/>simulation and target data
    
    loop Until fit is sufficient
        FiberPy->>+FiberCpp: 
        FiberCpp->>+FiberPy: 
        note right of FiberPy: FiberPy runs simulation,<br/>evaluates fit,<br/>saves progress,<br/>and adjusts model<br/>for next trial 
    end
-->

<img src = "images/fitting_sequence.png">

