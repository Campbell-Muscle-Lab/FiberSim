# -*- coding: utf-8 -*-
"""
Created on Thu May 20 16:54:41 2021

@author: kscamp3
"""

import numpy as np

def round_up_decade(val):
    # Returns the decade above val

    if (val > 0):
        y = np.power(10, (np.ceil(np.log10(val))))
    else:
        y = np.power(10, (np.ceil(np.log10(np.abs(val)))))
        
    return y

def round_down_decade(val):
    # Returns the decade below val
    
    if (val > 0):
        y = np.power(10, (np.floor(np.log10(val))))
    else:
        y = np.power(10, (np.floor(np.log10(np.abs(val)))))
        
    return y

def multiple_greater_than(val, mult):
    # Returns a multiple greater than

    return (mult * np.ceil(val/mult))

def multiple_less_than(val, mult):
    # Returns a multiple less than

    return (mult * np.ceil(val/mult))


def tidy_limits(vals, f=0.2):
    # Returns limits
    
    max_val = np.amax(vals)
    min_val = np.amin(vals)
    
    lims = np.NaN * np.ones(2)
    
    try:
        for (i, v) in enumerate([min_val, max_val]):
            if (i == 0):
                y = multiple_less_than(v, f*round_down_decade(v))
                if (y > 0):
                    if (y < 1):
                        d = int(np.abs(np.ceil(np.log10(y))-2))
                    else:
                        d=0
                else:
                    if (y > -1):
                        d = int(np.abs(np.ceil(np.log10(np.abs(y)))-2))
                    else:
                        d = 0
            else:
                y = multiple_greater_than(v, f*round_up_decade(v))
                if (y > 0):
                    if (y < 1):
                        d = int(np.abs(np.ceil(np.log10(y))-2))
                    else:
                        d=0
                else:
                    if (y > -1):
                        d = int(np.abs(np.ceil(np.log10(np.abs(y)))-2))
                    else:
                        d = 0
    
            y = np.around(y, decimals=d)
            
            lims[i] = y
    except:
        lims[0] = 0
        lims[1] = 1
      
    if ((lims[0] > 0) and (lims[1] > 0)):
        if (lims[0] < (0.7 * lims[1])):
            lims[0] = 0
        
    if ((lims[0] == lims[1])):
        lims[0] = lims[0] - 1
        lims[1] = lims[1] + 1

    return lims