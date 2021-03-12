---
title: Single run
has_children: false
parent: Getting started
grand_parent: Demos
nav_order: 1
---

# Single run

## Overview

This demo shows you how to run a simple simulation in FiberSim. You need a computer running Windows.

## What this demo does

This demo simulates a single half-sarcomere that is held isometric and activated in pCa 4.5 solution. Summary results are saved to a data file and plotted as a figure.

## Instructions

### Pre-requisites

If this is the first time using FiberSim, you need to:
+ [install software](../../installation/installation.html)

### Getting ready

Once you have installed FiberSim
+ Start Anaconda Navigator
+ Select the Environments tab (left-hand side)
+ Open a FiberSim terminal (left-click on the arrow-head to the right of FiberSim)
+ Change directory to `<repo>/code/FiberPy/FiberPy`, where `<repo>` is the directory you installed the software (e.g. `c:\temp\FiberSim`)

### Run a simulation

+ Type `python FiberPy.py "../../../demo_files/getting_started/single_run/batch_single_run.json"`
+ You should see some text appearing in the terminal window

## Viewing the results

+ Use Windows File Explorer to open `<repo>/demo_files/getting_started/single_run/simulation_output`
+ You should see
  + results.txt
  + summary.png
