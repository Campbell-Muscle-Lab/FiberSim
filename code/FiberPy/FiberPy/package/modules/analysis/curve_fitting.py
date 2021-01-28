# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:20:25 2020

@author: kscamp3
"""

import numpy as np
    
from scipy.optimize import curve_fit


def fit_pCa_data(x,y):
    """ Fits Hill-curve to x-y data """

    def y_pCa(x_data, pCa_50, n_H, y_min, y_amp):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = y_min + \
                y_amp * (np.power(np.power(10, -x), n_H) /
                (np.power(np.power(10, -x), n_H) + 
                     np.power(np.power(10, -pCa_50), n_H)))
        return y

    popt, pcov = curve_fit(y_pCa, x, y, [6.0, 2, np.amin(y), np.amax(y)])
    
    d = dict()
    d['pCa_50'] = popt[0]
    d['n_H']    = popt[1]
    d['y_min']  = popt[2]
    d['y_max']  = popt[3]
    d['x_fit']  = np.linspace(9.0, 4, 1000)
    d['y_fit']  = y_pCa(d['x_fit'], *popt)

    return d

def fit_hyperbola(x, y):
    """ Fits hyperbola of form (x+a)(y+b) = b*(x_0+a) to y data """
    
    def y_hyperbola(x_data, x_0, a, b):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = ((x_0+a)*b)/(x+a) - b
        return y
    
    popt, pcov = curve_fit(y_hyperbola, x, y,
                           [np.amax(x), 0.2*np.amax(x), 0.1])
    
    d = dict()
    d['x_0'] = popt[0]
    d['a'] = popt[1]
    d['b'] = popt[2]
    d['x_fit'] = np.linspace(0, np.amax(x), 1000)
    d['y_fit'] = y_hyperbola(d['x_fit'], *popt)
    
    return d

def fit_power_curve(x, y):
    """ Fits power curve of form y = x*b*(((x_0+a)/(x+a))-1) to y data """
    
    def y_power(x_data, x_0, a, b):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = x*b*(((x_0+a)/(x+a))-1)
        return y
    
    popt, pcov = curve_fit(y_power, x, y,
                           [np.amax(x), 0.2*np.amax(x), 0.1])
    
    d = dict()
    d['x_0'] = popt[0]
    d['a'] = popt[1]
    d['b'] = popt[2]
    d['x_fit'] = np.linspace(0, np.amax(x), 1000)
    d['y_fit'] = y_power(d['x_fit'], *popt)
    
    return d

def fit_exponential_recovery(x, y, n=1):
    """ Fits exponential recovery with a single exponential of form y = offset + amp*(1 - exp(-k*x)) to y data """
    
    if n==1:
        def y_single_exp(x_data, offset, amp, k):
            y = np.zeros(len(x_data))
            for i,x in enumerate(x_data):
                y[i] = offset + amp*(1 - np.exp(-k*x))
            return y
        
        popt, pcov = curve_fit(y_single_exp, x, y,
                               [y[0], y[-1]-y[0], (1/(0.2*np.amax(x)))])
        
        d = dict()
        d['offset'] = popt[0]
        d['amp'] = popt[1]
        d['k'] = popt[2]
        d['x_fit'] = np.linspace(x[0], x[-1], 1000)
        d['y_fit'] = y_single_exp(d['x_fit'], *popt)
        
        return d
        