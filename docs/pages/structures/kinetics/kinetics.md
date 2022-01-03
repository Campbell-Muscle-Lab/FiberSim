---
title: Kinetics
parent: Structures
has_children: false
nav_order: 7
---

# Kinetics

## Overview

The kinetic schemes for myosin and myosin-binding protein-C simulated by FiberCpp are defined in the [model file](../model/model.html).

Here is an example of a myosin scheme from the [isometric_activation_demo](../../demos/getting_started/isometric_activation/isometric_activation.html)

````
{
  "m_kinetics": [
    {
      "no_of_states": 3,
      "max_no_of_transitions": 2,
      "scheme": [
        {
          "number": 1,
          "type": "S",
          "extension": 0,
          "transition": [
            {
              "new_state": 2,
              "rate_type": "force_dependent",
              "rate_parameters": [
                10,
                200
              ]
            }
          ]
        },
        {
          "number": 2,
          "type": "D",
          "extension": 0,
          "transition": [
            {
              "new_state": 1,
              "rate_type": "constant",
              "rate_parameters": [
                200
              ]
            },
            {
              "new_state": 3,
              "rate_type": "gaussian",
              "rate_parameters": [
                82.5
              ]
            }
          ]
        },
        {
          "number": 3,
          "type": "A",
          "extension": 5.0,
          "transition": [
            {
              "new_state": 2,
              "rate_type": "poly",
              "rate_parameters": [
                75,
                1,
                2
              ]
            }
          ]
        }
      ]
    }
  ]
}
````

## Details

### Muscle

| Key | Comment |
| ---- | ---- |
| no_of_states | the number of states each myosin head can transitoin between |
| max_no_of_transitions | the maximum number of different transitions from a state in the model |

### Scheme

| Key | Comment |
| ---- | ---- |
| number | the state's number, typically 1 upwards |
| type | `S` super-relaxed, `D` disordered-relaxed, `A` attached |
| extension | the power-stroke distance for the state |

Super-relaxed heads can only exist as dimers. Thus, if a head transitions from an `S` to a `D` state, the partner head also transitons to `D`. Similarly, both heads for a dimer have to be in the `D` state in order for the pair to transition synchronusly to an `S` state.

### Transition

| Key | Comment |
| ---- | ---- |
| new_state | the new state number if the transition occurs |
| type | the type of rate, see below |
| rate parameters | variables defining the rate |

### Transition types

+ `constant`
  + rate_parameters:
    + a, s<sup>-1</sup> 
  + rate = a

+ `gaussian`
  + rate_parameters:
    + a, s<sup>-1</sup> 
  + rate = a e<sup>-k x<sup>2</sup></sup>

+ `force_dependent`
  + rate_parameters:
    + a, s<sup>-1</sup> 
    + b, s<sup>-1</sup> N<sup>-1</sup> m<sup>2</sup>
  + rate = a + b * node_force 

+ `poly`
  + rate_parameters
    + a, s<sup>-1</sup> 
    + b, depends on c
    + c, integer
	+ d, nm (optional)	
  + rate = a + b * (x + d)<sup>c</sup>
  
Note: if d is not specified, the code sets d = the state extension (5.0 in the myosin scheme example shown above)

+ `poly_asym`
  + rate_parameters
    + a, s<sup>-1</sup> 
    + b, depends on d
    + c, depends on d
	+ d, integer
	+ e, integer
	+ f, nm (optional)
  + rate = a + b * (x + f)<sup>d</sup> if x > f
  + rate = a + c * (x + f)<sup>e</sup> if x < f
  
Note: if f is not specified, the code sets f = the state extension (5.0 in the myosin scheme example shown above)
  


