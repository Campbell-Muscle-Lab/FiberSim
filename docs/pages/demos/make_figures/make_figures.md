---
title: Make figures
has_children: false
parent: Demos
nav_order: 4
---

# Make figures

## Overview

This demo shows how to make new figures from an existing results file.

## What this demo does

This demo uses the [output handler system](../../structures/output_handler.html) to make two figures from a [summary results file](../../structures/summary_results.html). You can build on this approach to make figures for presentations and papers.

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

+ Type `python FiberPy.py make_figures "../../../demo_files/getting_started/make_figures/output_handler.json"`
+ You should see some text appearing in the terminal window

## Viewing the results

+ Use Windows File Explorer to open `<repo>/demo_files/getting_started/make_figures/figures`
+ You should see
  + summary.png
  + sarah.png
