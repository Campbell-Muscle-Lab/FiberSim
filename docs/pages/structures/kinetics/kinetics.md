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

#### Constant

+ `constant`
  + rate_parameters:
    + a, s<sup>-1</sup> 
  + rate = a

#### Gaussian

+ `gaussian`
  + rate_parameters:
    + a, s<sup>-1</sup> 
  + rate = a e<sup>(-0.5 * m_k_cb * x<sup>2</sup>/(1e18 * k<sub>B</sub> * T))</sup>
  
  where m_k_cb is the crossbridge stiffness and T is the temperature (both defined in the [model file](../model/model.html)). k<sub>B</sub> is the Boltzmann constant. 
  
#### Gaussian_pc

  + `gaussian_pc`
  + rate_parameters:
    + a, s<sup>-1</sup> 
    + rate = a e<sup>(-0.5 * c_k_stiff * x<sup>2</sup>/(1e18 * k<sub>B</sub> * T))</sup>
  
  where c_k_stiff is the MyBP-C stiffness and T is the temperature (both defined in the [model file](../model.html)). k<sub>B</sub> is the Boltzmann constant. 

#### Force-dependent

+ `force_dependent`
  + rate_parameters:
    + a, s<sup>-1</sup> 
    + b, s<sup>-1</sup> N<sup>-1</sup> m<sup>2</sup>
  + rate = a + b * node_force 

#### Poly

+ `poly`
  + rate_parameters
    + a, s<sup>-1</sup> 
    + b, depends on c
    + c, integer
	+ d, nm (optional)	
  + rate = a + b * (x + d)<sup>c</sup>
  
Note: if d is not specified, the code sets d = the state extension (5.0 in the myosin scheme example shown above)

#### Poly_asym

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

#### Exp

+ `exp`
  + rate_parameters
    + a, s<sup>-1</sup> 
    + b, s<sup>-1</sup>
    + c, nm<sup>-1</sup> 
	+ d, nm
	+ e, nm
  + rate = a + b * exp<sup>(-c * (x + d))</sup> if x < e
  + rate = max_rate if x > e
  
  max_rate is defined in the [options file](../options/options.html).
  
#### Exp_wall

+ `exp_wall`
  + rate_parameters
    + a, s<sup>-1</sup> 
    + b, nm
    + c, nm
	+ d, nm<sup>-1</sup> 

  + rate = a * exp(- b * m_k_cb * (x + x_ext)/(1e18 * k<sub>B</sub> * T)) + max_rate * 1/(1 + exp(- d * (x - c)) 
  
  where m_k_cb is the crossbridge stiffness and T is the temperature (both defined in the [model file](../model/model.html)). k<sub>B</sub> is the Boltzmann constant.   max_rate is defined in the [options file](../options/options.html).

#### Exp_wall_sweep

+ `exp_wall_sweep`
  + rate_parameters
    + a, s<sup>-1</sup> 
    + b, nm
    + c, nm
	+ d, nm<sup>-1</sup>
    + e, s<sup>-1</sup> 

  + rate = a * exp(- b * m_k_cb * (x + x_ext)/(1e18 * k<sub>B</sub> * T)) + max_rate * 1/(1 + exp(- d * (x - c)) + 
  e * (2.0 - `no_of_active_neighbors`) 
  
  where m_k_cb is the crossbridge stiffness and T is the temperature (both defined in the [model file](../model/model.html)), k<sub>B</sub> is the Boltzmann constant, max_rate is defined in the [options file](../options/options.html).
  
  `no_of_active_neighbors` is the number of adjacent regulatory units on thin filament strand which are active. This number can be 0, 1, or 2, with the lower numbers being more likely when the muscle is minimally activated. If the sweep parameter `e` is greater than 0, myosin heads are more likely to detach when adjacent binding sites are inactive. This mimics the tropomyosin molecule 'sweeping' myosin heads off their binding sites.
  



