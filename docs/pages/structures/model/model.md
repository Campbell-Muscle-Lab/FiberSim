---
layout: default
title: Model
has_children: false
parent: Structures
nav_order: 1
---

# Model

In FiberSim, a `model file` defines the structure and the biophysical proporties of a the sarcomeric system that will be simulated.

| Level | Parameter | Units | Comments |
| ----  | ----      | ----  | ----     |
| muscle | no_of_half_sarcomeres | None | Integer defining the number of half-sarcomeres in series |
|  | no_of_myofibrils | None | Currently, must be 1 |
|  | sc_k_stiff | N m<sup>-1</sup> | Stiffness of an optional series elastic component. Assumed to be infinite if not defined |
|  | initial_hs_length | nm | Initial length of half-sarcomere(s) |
|  | prop_fibrosis | None | Proportion of cross-section occupied by fibrosis |
|  | prop_myofilaments | None | Proportion of cross-section that is not fibrosis that is occupied by myofilaments |
|  | m_filament_density | m<sup>-2</sup> | The number of thick filaments per m<sup>2</sup> in a myofilament |


| Level | Parameter | Units | Comments |
| ----  | ----      | ----  | ----     |
| lattice_parameters | viscosity | N m<sup>-2</sup> (nm half-sarcomere<sup>-1</sup>)<sup>-1</sup> | Passive viscous drag between myofilaments |

| Level | Parameter | Units | Comments |
| ----  | ----      | ----  | ----    |
| thick_structure | m_n | None | Integer defining the number of thick filaments to simulate. Must be an integer squared (4, 9, 16, etc. up to 196). Larger numbers produce smoother traces but take longer to simulate.|
|  | m_crowns_per_filament | None | Integer defining the number of crowns per thick filament |
|  | m_myosin_hubs_per_crown | None | Integer defining the number of hubs per crown. Currently must be 2 for myosin dimers. |
|  | m_inter_crown_rest_length | nm | Distance between crowns in an unstressed thick filament |
|  | m_lambda | nm | Length of the thick filament bare zone |
|  | m_starting_angle | degrees | Orientation of the first hub |
|  | m_inter_crown_twist | degrees | Angular twist between crowns |
|  | m_within_hub_twist | degrees | Angular twist between heads in a hub |


