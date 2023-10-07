---
layout: default
title: Demos
has_children: True
nav_order: 10
---

# Demos
{:.no_toc}

This section of the documentation demonstrates how to use FiberSim.

A high-level overview is provided immediately below. The sub-pages provide examples that show how to:

+ simulate single trials
+ calculate force-pCa curves with optional step-length changes and/or k<sub>tr</sub> maneuvers
+ calculate force-velocity and force-power curves
+ compare models
+ simulate myofibrils composed of half-sarcomeres connected in series
+ run simulations with parameters selected from a defined parameter space

## Note

None of the simulations in the demonstrations have been fitted to data. They are designed to help people use the software and not intended to reproduce specific experimental results.

# Overview

The overarching strategy is to automate as much as possible.

As shown below, users can run every demo on this site by launching FiberPy with the following files:
+ a model file
+ a setup file
+ an options file

FiberPy handles everything else.

The model file describes the properties of a *base* half-sarcomere and, in the case of myofibrils, how many of them are arranged in series.

The setup file describes the experimental protocol (e.g. pCa levels, length-changes, loadeded shortening) and provides options for comparing simulations at different starting lengths, models with different parameters, etc.

The options file provides fine control over the way FiberCpp runs simulations and, if desired, how FiberCpp will write status files to disk.

The demos explain how to use most of FiberSim's current capabilities. If you want to try and use FiberSim to do something that is not documented, submit a request or contact us.

<img src="images/FiberSim_workflow.png">