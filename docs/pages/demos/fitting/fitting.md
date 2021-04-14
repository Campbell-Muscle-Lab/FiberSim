---
title: Fitting
has_children: True
parent: Demos
nav_order: 7
---

## Fitting

In this context, fitting means adjusting model parameters to obtain the best possible fit between the simulations and experimental data.

This can be done by fitting:
+ in the time domain
  + matching simulations and experimental traces plotted against time
+ pCa curves
  + matching simulations and experimental data plotted as a function of activating Ca<sup>2+</sup> concentration
+ in the frequency domain
  + matching Nyquist plots showing elastic and viscous moduli as a function of frequency


### General strategy

The general strategy in pseudo-code is as follows:

````
define a starting_model
set the working_model equal to the starting_model
while (fit_is_not_good_enough)
{    
    run a simulation using the working_model
    evaluate fit by comparing the simulation to a target 
    adjust the working_model by tweaking the parameters
}
````

### Optimization structure

Fitting in MATMyoSim is performed using an [Optimization structure](..\..\structures\optimization_structure\optimization_structure.html)

