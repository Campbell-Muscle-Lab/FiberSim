---
layout: default
title: Code Validation
has_children: true
nav_order: 6
---

# Code Validation
{:.no_toc}

At each time step, FiberCpp calculates the position of each thick and thin filaments node. FiberCpp also computes the kinetics for the actin binding sites and the myosin/myosin-binding protein C molecules. In order to make sure that the calculations are done properly, we have developed a testing suite in Python.

## Force-balance Test

As explained [here](../FiberCpp/calculations/calculations.html), FiberCpp solves a matrix equation of form $K x = F$, where $x$ is a vector containing the node positions, $K$ is the tridiagonal stiffness matrix, and $F$ is a vector containing the cross-bridges, titin, and myosin-binding protein C forces. The force-balance test is a Python code written to evaluate the accuracy of the algorithm calculations.

## Kinetics Test

FiberCpp also computes the kinetics for the actin binding sites, the myosin and the myosin-binding protein C molecules. At each time step, transition probabilities are calculated according to the rate laws implemented in the [model file](../structures/model/model.html) and transition events are implemented. 

The kinetics test is a Python code written to make sure that the transition events calculated by FiberCpp occur according to the rate laws that were specified in the model file. The kinetics test suite validates actin, myosin, and myosin-binding protein C (MyBPC) kinetics.
