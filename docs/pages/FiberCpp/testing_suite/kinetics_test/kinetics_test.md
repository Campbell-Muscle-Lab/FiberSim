---
title: Kinetics test
parent: Testing suite
grand_parent: FiberCpp
has_children: true
nav_order: 2
---

# Kinetics test 

FiberCpp also computes the kinetics for the actin binding sites, the myosin and the myosin-binding protein C molecules. At each time step, transition probabilities are calculated according to the rate laws implemented in the [model file](../../../structures/model/model.html) and transition events are implemented. 

The kinetics test is a Python code written to make sure that the transition events calculated by FiberCpp occur according to the rate laws that were specified in the model file. The kinetics test suite validates the thin and thick kinetics. 

## Thin filament kinetics 

Actin regulatory units are composed of binding sites that can transition from an inactive state to an active state, and vice-versa. At each time step, all the binding sites from a regulatory unit can transition together according to the rate constants `a_k_on` and `a_k_off` defined in the [model file](../../../structures/model/model.html). The probabilty of a regulatory unit activating is calculated using +a_{k_{on}} [Ca^{2+}]+, and the probability of a regulatory unit deactivating is calculated using +a_{k_{off}}+. These probabilities are also modulated by `a_k_coop`, which is a parameter quantifiying the inter-regulatory unit cooperativity:

+ If no direct neighboring unit is activated, the activation rate is given by: 

++ a_{k_{on}} [Ca^{2+}] ++  

+ If one direct neighboring unit is already activated, the activation rate is given by:

++ a_{k_{on}} [Ca^{2+}] (1 + a_{k_{coop}})++  

+ If the two direct neighboring units are already activated, the activation rate is given by:

++ a_{k_{on}} [Ca^{2+}] (1 + 2 \, a_{k_{coop}})++   

The same reasoning goes for the inactivating rate constant +a_{k_{off}}+.


### What this test does

The actin kinetics test:

+ Runs a simulations in which a half-sarcomere is held isometric and activated in a solution with a pCa of 6.0. 

+ Saves a status file at each time step, which contains the state of every regulatory unit. 

+ Assesses all the regulatory units state transitions occuring between two consecutive time-steps, and calculates the apparent rate constants.

+ Compares the calculated rate constants with those provided in the model file. 

The batch file containing the instructions for this test is:

```
{
  "FiberSim_batch": {
    "FiberCpp_exe": {
      "relative_to": "this_file",
      "exe_file": "../../bin/FiberCpp.exe"
    },
    "job": [
      {
        "relative_to": "this_file",
        "model_file": "sim_input/model.json",
        "options_file": "sim_input/options.json",
        "protocol_file": "sim_input/pCa6_protocol.txt",
        "results_file": "sim_output/results.txt",
        "output_handler_file": "sim_input/output_handler.json"
      }
    ],
    "batch_validation": [
      {
        "validation_type": "a_kinetics",
        "relative_to": "this_file",
        "model_file": "sim_input/model.json",
        "options_file": "sim_input/options.json",
        "protocol_file": "sim_input/pCa6_protocol.txt",
        "output_data_folder": "sim_output/analysis"
      }
    ]
  }
}
```

As explained above, this batch file consists in a single job. The block called `batch_validation` provides the type of validation that must be run and the necessary files information.


### Instructions

Before proceeding, make sure that you have followed the [installation instructions](../../../installation/installation.html) and that you already tried to run the [getting started demos](../../../demos/getting_started/getting_started.html).

### Getting ready

+ Open an Anaconda Prompt

+ Activate the FiberSim Anaconda Environment by executing:
```
conda activate fibersim
```
+ Change directory to `<FiberSim_dir>/code/FiberPy/FiberPy`, where `<FiberSim_dir>` is the directory where you installed FiberSim. 

### Run the test

+ Type:
```
python Fiberpy.py run_batch "../../../testing_suite/thin_kinetics/batch_a_kinetics.json"
```

+ You should see text appearing in the terminal window, showing that the simulations are running. When it finishes (this can take ~30 minutes), you should see something similar to the image below.

![command prompt](figures/command_prompt_a.png)

### Viewing the results

The results and summary figure from the simulation are written to files in `<FiberSim_dir>/testing_suite/thin_kinetics/sim_output`

<img src='figures/sim_output_thin.PNG'>

The `hs` folder contains the status files that were dumped at each time-step calculation.

The `analysis` folder contains two figure files showing the calculated activating and deactivating rate constants for the three cooperative scenarios, as well as the model values.

<img src='figures/active_to_inactive.png' width="48%">
<img src='figures/inactive_to_active.png' width="48%">

## Thick filament kinetics 

Thick filaments are composed of myosin and myosin-binding protein C. The kinetics scheme for both molecules is defined in the model file, and examples will be provided below.
The transition parameters between each states are described in the `transition` blocks. The two following tests respectively verifiy the myosin and the myosin-binding protein C kinetics. 

## Myosin kinetics test

### What this test does

The myosin kinetics test:

+ Runs a simulations in which a half-sarcomere is held isometric and activated in a solution with a pCa of 4.5. 

+ Saves a status file at each time step, which contains the state of every myosin molecules. 

+ Assesses all the myosin transitions occuring between two consecutive time-steps, and calculates the apparent rate constants.

+ Compares the calculated rate constants with those provided in the following model file:

```
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
                100,
                100
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
                200
              ]
            }
          ]
        },
        {
          "number": 3,
          "type": "A",
          "extension": 7.0,
          "transition": [
            {
              "new_state": 2,
              "rate_type": "poly",
              "rate_parameters": [
                150,
                1,
                2
              ]
            }
          ]
        }
      ]
    }
  ]
```


### Getting ready

+ Open an Anaconda Prompt

+ Activate the FiberSim Anaconda Environment by executing:
```
conda activate fibersim
```
+ Change directory to `<FiberSim_dir>/code/FiberPy/FiberPy`, where `<FiberSim_dir>` is the directory where you installed FiberSim. 

### Run the test

+ Type:
 ```
python Fiberpy.py run_batch "../../../testing_suite/thick_kinetics/myosin/batch_m_kinetics.json"
 ```

+ You should see text appearing in the terminal window, showing that the simulations are running. When it finishes (this can take ~20-30 minutes), you should see something similar to the image below.

![command prompt](figures/command_prompt_m.PNG)

### Viewing the results

The results and summary figure from the simulation are written to files in `<FiberSim_dir>/testing_suite/thick_kinetics/myosin/sim_output`

<img src='figures/sim_output_myosin.PNG'>

The `hs` folder contains the status files that were dumped at each time-step calculation.

The `analysis` folder contains figure files showing the calculated rate laws for all myosin transitions, as well as the model rate laws.

Transition #0 (from state 1 to state 2) is force-dependent, and the transition rate is plotted as a function of the node force.

<img src='figures/iso_0_transition_0.png' width="50%">

Transition #1 (from state 2 to state 1) is constant.

<img src='figures/iso_0_transition_1.png' width="50%">

Transition #2 (from state 2 to state 3) is stretch-dependent, and the transition rate is plotted as a function of the stretch values. Furthermore, the attachement probabilities are modulated by the angle difference between a myosin head and the nearest binding sites. Each panel represents an angle difference interval and shows the calculated and theoretical attachment rates for the mid-interval angle difference. 

<img src='figures/iso_0_transition_2.png' width="100%">

Transition #3 (from state 3 to state 2) is stretch-dependent, and the transition rate is plotted as a function of the stretch values.

<img src='figures/iso_0_transition_3.png' width="50%">

## Myosin-binding protein C kinetics test

### What this test does

The mybpc kinetics test:

+ Runs a simulations in which a half-sarcomere is held isometric and activated in a solution with a pCa of 4.5. 

+ Saves a status file at each time step, which contains the state of every mybpc molecule. 

+ Assesses all the mybpc transitions occuring between two consecutive time-steps, and calculates the apparent rate constants.

+ Compares the calculated rate constants with those provided in the following model file:

```
 "c_kinetics": [
    {
      "no_of_states": 2,
      "max_no_of_transitions": 1,
      "scheme": [
        {
          "number": 1,
          "type": "D",
          "extension": 0,
          "transition": [
            {
              "new_state": 2,
              "rate_type": "gaussian_pc",
              "rate_parameters": [
                200
              ]
            }
          ]
        },
        {
          "number": 2,
          "type": "A",
          "extension": 0.0,
          "transition": [
            {
              "new_state": 1,
              "rate_type": "poly",
              "rate_parameters": [
                150,
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
```

### Getting ready

+ Open an Anaconda Prompt

+ Activate the FiberSim Anaconda Environment by executing:
```
conda activate fibersim
```
+ Change directory to `<FiberSim_dir>/code/FiberPy/FiberPy`, where `<FiberSim_dir>` is the directory where you installed FiberSim. 

### Run the test

+ Type:
 ```
python Fiberpy.py run_batch "../../../testing_suite/thick_kinetics/mybpc/batch_c_kinetics.json"
 ```

+ You should see text appearing in the terminal window, showing that the simulations are running. When it finishes (this can take ~20-30 minutes), you should see something similar to the image below.

![command prompt](figures/command_prompt_c.PNG)

### Viewing the results

The results and summary figure from the simulation are written to files in `<FiberSim_dir>/testing_suite/thick_kinetics/mybpc/sim_output`

<img src='figures/sim_output_mybpc.PNG'>

The `hs` folder contains the status files that were dumped at each time-step calculation.

The `analysis` folder contains figure files showing the calculated rate laws for all mybpc transitions, as well as the model rate laws.

<img src='figures/iso_0_transition_0_mybpc.png' width="100%">
<img src='figures/iso_0_transition_1_mybpc.png' width="50%">






