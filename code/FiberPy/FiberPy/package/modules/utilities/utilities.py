# -*- coding: utf-8 -*-
"""
Created on Thu May 20 16:54:41 2021

@author: kscamp3
"""

import numpy as np

def round_up_decade(val):
    # Returns the decade above val

    return np.power(10,(np.ceil(np.log10(val))))

def multiple_greater_than(val, mult):
    # Returns a multiple greater than

    return (mult * np.ceil(val/mult))

def tidy_limits(vals, f=0.2):
    # Returns limits

    max_val = np.amax(vals)
    
    top = multiple_greater_than(max_val, f*round_up_decade(max_val))
    if (top < 1):
        d = int(np.abs(np.ceil(np.log10(top))-2))
    else:
        d=0
    top = np.around(top, decimals=d)

    return [0, top]