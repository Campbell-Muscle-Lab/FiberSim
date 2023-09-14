# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:20:25 2020

@author: kscamp3
"""

import numpy as np
    
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

import matplotlib.pyplot as plt


def fit_pCa_data(x,y):
    """ Fits Hill-curve to x-y data """

    def y_pCa(x_data, pCa_50, n_H, y_min, y_amp):
        
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
           
            x = float(x)
            y[i] = y_min + \
                y_amp * (np.power(np.power(10, -x), n_H) /
                (np.power(np.power(10, -x), n_H) + 
                     np.power(np.power(10, -pCa_50), n_H)))
        return y

    try:
        min_bounds = [3.0, 0.01, -np.inf, 0]
        max_bounds = [10.0, 100.0, np.inf, np.inf]
        popt, pcov = curve_fit(y_pCa, x, y,
                               [6.0, 2, np.amax([0, np.amin(y)]), np.amax([0, np.amax(y)])],
                               bounds=(min_bounds, max_bounds))
    except:
        print('fit_pCa_data failed')
        popt = [6.1, 2.1, 0.1, 1.0]

    d = dict()
    d['pCa_50'] = popt[0]
    d['n_H']    = popt[1]
    d['y_min']  = popt[2]
    d['y_amp']  = popt[3]
    d['x_fit']  = np.linspace(9.0, 4, 1000)
    d['y_fit']  = y_pCa(d['x_fit'], *popt)
    d['y_predict'] = y_pCa(x, *popt)

    return d

def fit_IC_50(x,y, type_curve = "increasing"):
    """ Fits increasing or decreasing Hill-curve to x-y data for drug-response"""

    def y_drug_pos(x_data, IC_50, n_H, y_min, y_amp): # increasing Hill curve
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = y_min + \
                y_amp * (np.power(x, n_H) /
                (np.power(x, n_H) + 
                     np.power(IC_50, n_H)))
        return y

    def y_drug_neg(x_data, IC_50, n_H, y_min, y_amp): # decreasing Hill curve
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = (y_amp + y_min) - \
                y_amp * (np.power(x, n_H) /
                (np.power(x, n_H) + 
                     np.power(IC_50, n_H)))
        return y

    try:
        
        if type_curve == "increasing":
            popt, pcov = curve_fit(y_drug_pos, x, y,
                    [0.5, 1.5, np.amax([0, np.amin(y)]), np.amax([0, np.amax(y)])])
        elif type_curve == "decreasing":
            popt, pcov = curve_fit(y_drug_neg, x, y,
                    [0.5, 1.5, np.amax([0, np.amin(y)]), np.amax([0, np.amax(y)])])            
    except:
        print('fit_IC_50 failed')
        popt = [0.5, 1.5, np.amax([0, np.amin(y)]), np.amax([0, np.amax(y)])]                               

    d = dict()
    d['IC_50'] = popt[0]
    d['n_H']    = popt[1]
    d['y_min']  = popt[2]
    d['y_max']  = popt[3]
    d['x_fit']  = np.linspace(0.01, 100, 10000)
    if type_curve == "increasing":
        d['y_fit']  = y_drug_pos(d['x_fit'], *popt)
    elif type_curve == "decreasing":
        d['y_fit']  = y_drug_neg(d['x_fit'], *popt)

    return d

def fit_hyperbola(x, y):
    """ Fits hyperbola of form (x+a)(y+b) = b*(x_0+a) to y data """
    
    def y_hyperbola(x_data, x_0, a, b):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = ((x_0+a)*b)/(x+a) - b
        return y
    
    try:
        popt, pcov = curve_fit(y_hyperbola, x, y,
                           [np.amax(x), 0.2*np.amax(x), 0.3])

    except:

        print('fit hyperbola failed')
        popt = [np.amax(x), 0.2*np.amax(x), 0.3]  
    
    d = dict()
    d['x_0'] = popt[0]
    d['a'] = popt[1]
    d['b'] = popt[2]
    d['x_fit'] = np.linspace(0, np.amax(x), 1000)
    d['y_fit'] = y_hyperbola(d['x_fit'], *popt)
    d['y_predict'] = y_hyperbola(x, *popt)
    
    return d

def fit_power_curve(x, y):
    """ Fits power curve of form y = x*b*(((x_0+a)/(x+a))-1) to y data """
    
    def y_power(x_data, x_0, a, b):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = x*b*(((x_0+a)/(x+a))-1)
        return y
    
    try:
        popt, pcov = curve_fit(y_power, x, y,
                           [np.amax(x), 0.2*np.amax(x), 0.1])
    except:
        print('fit_power_curve failed')
        popt = [np.amax(x), 0.2*np.amax(x), 0.1]
        
    
    d = dict()
    d['x_0'] = popt[0]
    d['a'] = popt[1]
    d['b'] = popt[2]
    d['x_fit'] = np.linspace(0, np.amax(x), 1000)
    d['y_fit'] = y_power(d['x_fit'], *popt)
    
    return d

def fit_exponential(x, y, n=1):
    
    # Test whether system is building or declining
    # lin_mod = fit_straight_line(x,y)
    if (y[-1] > y[0]):
        build_mode = 1
    else:
        build_mode = -1

    sorted_y = np.sort(y)
    min_y = np.amin(y)
    max_y = np.amax(y)
    
    if (build_mode == 1):
        guess_half_index = np.argmax(sorted_y > (min_y + (0.5 * (max_y - min_y))))
        p = [max_y, (min_y - max_y)]
    else:
        y_sorted = sorted_y[::-1]
        guess_half_index = len(y_sorted) - np.argmax(y_sorted < (min_y + 0.5 * (max_y - min_y))) - 1
        p = [min_y, (max_y-min_y)]
        
    k = -np.log(0.5) / (x[guess_half_index] - x[0])
    p.append(k)
    
    def y_func(x_data, offset, amp, k):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = offset + amp*(-abs(k) * x)
        return y
    
    min_bounds = [-np.inf, -np.inf, 0.0]
    max_bounds = [np.inf, np.inf, np.inf]
    
    popt, pcov = curve_fit(y_func, x, y,
                           p,
                           bounds=(min_bounds, max_bounds),
                           maxfev=5000)
    
    d = dict()
    d['offset'] = popt[0]
    d['amp'] = popt[1]
    d['k'] = popt[2]
    d['x_fit'] = x
    d['y_fit'] = y_func(d['x_fit'], *popt)
    d['r_squared'] = r2_score(y, d['y_fit'])
    
    return d
    
    
    

def fit_exponential_recovery(x, y, n=1):
    """ Fits exponential recovery with a single exponential of form y = offset + amp*(1 - exp(-k*x)) to y data """
    
    if n==1:
        def y_single_exp(x_data, offset, amp, k):
            y = np.zeros(len(x_data))
            for i,x in enumerate(x_data):
                y[i] = offset + amp*(1 - np.exp(-k*x))
            return y

        min_bounds = [-np.inf, -np.inf, 0.0]
        max_bounds = [np.inf, np.inf, np.inf]
        
        

        popt, pcov = curve_fit(y_single_exp, x, y,
                               [y[0], y[-1]-y[0], (1/(0.2*np.amax(x)))],
                               bounds=(min_bounds, max_bounds),
                               maxfev=5000)
        
        d = dict()
        d['offset'] = popt[0]
        d['amp'] = popt[1]
        d['k'] = popt[2]
        d['x_fit'] = x
        d['y_fit'] = y_single_exp(d['x_fit'], *popt)
        d['r_squared'] = r2_score(y, d['y_fit'])
        
        return d

def fit_exponential_decay(x, y):
    """ Fits exponential decay with a single exponential of form y = offset + amp*exp(-k*x) to y data """    

    def y_single_exp(x_data, offset, amp, k):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = offset + amp*np.exp(-k*(x))
        return y   

    min_bounds = [-np.inf, -np.inf, 0.0]
    max_bounds = [np.inf, np.inf, np.inf]
    
    st = [y[0], y[-1] - y[0], -np.log(0.5)/0.5*(np.amax(x)+np.amin(x))]
    st = [1091, 40.24, 0.3]
    print(st)
    
    print(x[0])
    
    try:
        
        popt, pcov, info, mesg = curve_fit(y_single_exp, x, y, [y[-1], y[0]-y[-1], -np.log(0.5)/0.5*(np.amax(x)+np.amin(x))],
                               maxfev=1000,
                               bounds=(min_bounds, max_bounds))
        
        print(info)
        print(mesg)
        
    except:
        print('fit exponential decay failed - setting decay rate to NaN')
        popt = [y[-1], y[0]-y[-1], np.nan]

    d = dict()
    d['offset'] = popt[0]
    d['amp'] = popt[1]
    d['k'] = popt[2]
    d['x_fit'] = np.linspace(x[0], x[-1], 1000)
    d['y_fit'] = y_single_exp(d['x_fit'], *popt)
        
    return d

def fit_straight_line(x, y):
    """ Fits a straight line to data """

    import statsmodels.api as sm

    # Create a regression model
    x = sm.add_constant(x)
    model  = sm.OLS(y, x)

    results = model.fit()

    # Create a dictionary for the results
    d = dict()
    d['x'] = x
    d['y'] = y
    d['intercept'] = results.params[0]
    d['slope'] = results.params[1]

    return d
