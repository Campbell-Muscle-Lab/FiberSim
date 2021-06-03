# -*- coding: utf-8 -*-
"""
Created on Thu May 20 16:54:41 2021

@author: kscamp3
"""

import numpy as np

def multiple_greater_than(val, mult):
    # Returns a multiple greater than

    return (mult * np.ceil(val/mult))