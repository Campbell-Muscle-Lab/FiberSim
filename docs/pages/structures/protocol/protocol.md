---
title: Protocol
parent: Structures
has_children: false
nav_order: 20
---

# Protocol

## Overview

The protocol file for a [job](../job/job.html) defines the experiment that is being simulated.

It is a tab-delimited text file with 4 columns and one row for each time-step in the simulation. Here is a very simple example.

````
dt	pCa	dhsl	mode
0.0001	9.0	0.0	-2.0
0.0001	9.0	0.0	-2.0
0.0001	4.5	0.0	-2.0
0.0001	4.5	0.0	-2.0
0.0001	4.5	1.0	-2.0
0.0001	4.5	0.0	-2.0
0.0001	4.5	0.0	1000.0
0.0001	4.5	0.0	1000.0
````

## Details

This section explains the components

| Key | Comment |
| ---- | ---- |
| dt | the duration of the time-step in seconds |
| pCa | the Ca<sup>2+</sup> concentration |
| dhsl | the length change in nm per half-sarcomere imposed during the time-step |
| mode | explained below |

## Mode

The `mode` value for each time-step controls the loading condition for the time-step.

+ -2 - the simulation is in length control mode
  + if `dhsl` is
    + 0, the simulation is isometric
    + >0, the muscle is extended
    + <0, the muscle is shortened

+ -1 - slack mode
  + this means check whether the muscle has fallen slack (that is, force < 0>)
    + if it is slack, shorten against zero load
    + if not, change length as defined by `dhsl`

+ >=0 - isotonic mode
  + lengthen or shorten the muscle as required to keep stress equal to the mode value

## Example

````
dt	pCa	dhsl	mode
0.0001	9.0	0.0	-2.0
0.0001	9.0	0.0	-2.0
0.0001	4.5	0.0	-2.0
0.0001	4.5	0.0	-2.0
0.0001	4.5	1.0	-2.0
0.0001	4.5	0.0	-2.0
0.0001	4.5	0.0	1000.0
0.0001	4.5	0.0	1000.0
````

This example instruction file (copied from the top of the page) defines a simulation of 0.8 ms (there are 8 time-steps, each with a `dt` of 0.1 ms).

The first two data lines show that the muscle is held isometric (`mode` is -2.0, `dhsl` is 0).

The muscle is then activated in pCa 4.5 solution for the remainder of the simulation (`pCa` is 4.5 for time-steps 3 onwards).

The muscle is stretched by 1 nm per half-sarcomere on step 5 (`dhsl` is 1) and held at an isotonic load of 1000 N m<sup>-2</sup> for the time-steps 7 and 8 (`mode` is 1000).

## More details

See the [demos](../../demos/demos.html)
