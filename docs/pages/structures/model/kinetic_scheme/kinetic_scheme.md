---
title: Kinetic scheme
has_children: false
parent: Model
grand_parent: Structures
nav_order: 1
---

# Kinetic scheme

Here is an example

````
"m_kinetics": {
        "no_of_states": 2,
        "max_no_of_transitions": 2,
        "scheme": [
        {
          "number": 1,
          "type": "D",
          "extension": 0,
          "transition":
          [
           {
            "new_state": 2,
            "rate_type": "gaussian",
            "rate_parameters": [100, 0, 3]
           }
          ]
        },
        {
          "number": 2,
          "type": "A",
          "extension": 7.0,
          "transition":
          [
           {
            "new_state": 1,
            "rate_type": "poly",
            "rate_parameters": [100, 5, 4]
           }
        ]
      }
    ]
  }
````
