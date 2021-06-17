---
title: Results
parent: Structures
has_children: false
nav_order: 25
---

# Results

## Overview

FiberCpp saves a summary of the simulation results. The file name is set in the [job](../job/job.html) structure.

The result file provides a high-level overview of the simulation and includes metrics such as force, half-sarcomere length, and the relative populations of different actin and myosin states.

Additional information about the status of individual molecules can be saved in [status files](../status_file/status_file.html).

## Details

Results files are saved as tab-delimited text files.

[Click here for an example file](example_result_file.txt)

Here are the top few lines

````
time	pCa	hs_length	hs_command_length	hs_slack_length	force	titin_force	extracellular_force	a_fil_length	m_fil_length	a_pop_0	a_pop_1	m_pop_0	m_pop_1	m_pop_2	c_pop_0	c_pop_1
0.0001	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.995756	0.00424383	0	1	0
0.0002	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.994599	0.00540123	0	1	0
0.0003	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.991512	0.00848765	0	1	0
0.0004	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.989583	0.0104167	0	1	0
0.0005	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.987269	0.0127315	0	1	0
0.0006	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.984568	0.0154321	0	1	0
0.0007	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.98341	0.0165895	0	1	0
0.0008	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.982639	0.0173611	0	1	0
0.0009	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.98071	0.0192901	0	1	0
0.001	9.000	1000	1000	-nan(ind)	2861.55	2861.55	0	1015.88	809.004	1	0	0.978781	0.0212191	0	1	0
````

| Key | Comment |
| ---- | ---- |
| time | the time in s |
| pCa | the Ca<sup>2+</sup> concentration |
| hs_length | the length of the half-sarcomere in nm |
| command_length | the command length of the half-sarcomere in nm, will be equal to command length unless the muscle is slack |
| slack_length | the length of the half-sarcomere in nm when force is zero, shows -nan(ind) if not calculated for the time-step |
| force | force per unit area in N m<sup>-2</sup> |
| titin_force | titin-based force per unit area in N m<sup>-2</sup> |
| extracellular_force | extracellular_force per unit area in N m<sup>-2</sup> |
| a_fil_length | average length of thin filaments in nm |
| m_fil_length | average length of thick filaments in nm |
| a_pop_i | the relative population of binding sites in state i |
| m_pop_i | the relative population of heads in state i |
| c_pop_i | the relative population of myosin binding protein-C molecules in state i |
