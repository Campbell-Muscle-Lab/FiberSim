---
title: Single run
has_children: false
parent: Getting started
grand_parent: Demos
nav_order: 1
---

# Single run

## Overview

This demo shows you how to run a single simulation in FiberSim.

If you read the documentation, you will see that FiberSim is composed of several parts
+ the [core model](../core_model/core_model.html)
  + which implements the calculations underlying the simulation
  + is written in C++
  + is a console application stored as FiberSim.exe in `repo/bin`
+ the [interface](../inteface/interface.html)
  + which provides functions that simplify
    + simplify running simulations
    + analyzing the output
    + fitting models to experimental data
  + is written in Python

## Instructions

+ Start Anaconda Navigator
+ Select the Environments tab (left-hand side)
+ Open a FiberSim terminal
+ Change directory to `repo/code/python/FiberSim_utilities`
+ Type `python FiberSim_utilities.py demos getting_started single_run`
+ You should see
  + some text appearing in the terminal window
  + a new figure popping up
+ Close the figure to return focus to the terminal window

## Video

Click on the screenshot below for a video demo

<a href="https://drive.google.com/file/d/1IqP5XdBfmSc9TSxgKWQoyXUfLmr64CJ4/view?usp=sharing">
![Video screenshot](single_run_screenshot.png)</a>






