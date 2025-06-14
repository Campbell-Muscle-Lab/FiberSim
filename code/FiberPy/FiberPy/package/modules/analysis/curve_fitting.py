# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:20:25 2020

@author: kscamp3
"""

import numpy as np

from numpy.polynomial import Polynomial as poly
    
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar
from sklearn.metrics import r2_score

import matplotlib.pyplot as plt


def fit_pCa_data(x, y, no_of_fit_points = 1000):
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
        min_bounds = [4.0, 0.01, -np.inf, 0]
        max_bounds = [8.0, 100.0, np.inf, np.inf]
        popt, pcov = curve_fit(y_pCa, x, y,
                               [6.0, 2, np.amax([0, np.amin(y)]), np.amax([0, np.amax(y)])],
                               bounds=(min_bounds, max_bounds))
    except:
        print('fit_pCa_data failed')
        popt = [6.0, 10.0, 0, 1e5]

    d = dict()
    d['pCa_50'] = popt[0]
    d['n_H']    = popt[1]
    d['y_min']  = popt[2]
    d['y_amp']  = popt[3]
    d['x_fit']  = np.linspace(9.0, 4, no_of_fit_points)
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
    
    def neg_y_power(x, x_0, a, b):
        y = -y_power(np.asarray([x]), x_0, a, b)
        return y
    
    try:
        popt, pcov = curve_fit(y_power, x, y,
                           [np.amax(x), 0.2*np.amax(x), np.amax(y) / 0.1 * np.amax(x)])
    except:
        print('fit_power_curve failed')
        popt = [np.amax(x), 0.2*np.amax(x), 0.1]
        
    # Get max of curve
    r = minimize_scalar(neg_y_power, bounds=[0, np.amax(x)],
                       args=(popt[0],popt[1],popt[2]))
    
    d = dict()
    d['x_0'] = popt[0]
    d['a'] = popt[1]
    d['b'] = popt[2]
    d['x_at_max_power'] = r['x'][0]
    d['max_power'] = -r['fun'][0]
    d['x_fit'] = np.linspace(0, np.amax(x), 1000)
    d['y_fit'] = y_power(d['x_fit'], *popt)
    d['y_predict'] = y_power(x, *popt)
    
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
    
    def y_func(x_data, offset, amp, k2):
        y_guess = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y_guess[i] = offset + amp*(-abs(k2) * x)
        return y_guess
    
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
    
    d['ken'] = y_func(np.linspace(0, 1, 10), 1, 3,-0.1)

    print(d)
    print(*popt)
    exit(1)
    
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
        
        
        try:
            popt, pcov = curve_fit(y_single_exp, x, y,
                               [y[0], y[-1]-y[0], (1/(0.2*np.amax(x)))],
                               bounds=(min_bounds, max_bounds),
                               maxfev=5000)
        except:
            print('fit exponential decay failed - setting decay rate to bad values')
            popt = [1e5, 1e5, 10000]

        
        d = dict()
        d['offset'] = popt[0]
        d['amp'] = popt[1]
        d['k'] = popt[2]
        d['x_fit'] = x
        d['y_fit'] = y_single_exp(d['x_fit'], *popt)
        d['r_squared'] = r2_score(y, d['y_fit'])
        
        return d
    
def fit_shortening_length_trace(x, y):
    """ Fits shortening length trace """
    
    def y_function(x, a, b, c):
        y_data = np.zeros(len(x))
        for i, x_point in enumerate(x):
            y_data[i] = a + (b * np.exp(-c*x_point))
        return y_data
            
    s = [1000, 100, 10]
    min_bounds = [0, 0, 0]
    max_bounds = [np.inf, np.inf, 100]
    
    popt, pcov = curve_fit(y_function, x, y, s,
                           bounds = (min_bounds, max_bounds),
                           max_nfev=5000)
    
    d = dict()
    d['a'] = popt[0]
    d['b'] = popt[1]
    d['c'] = popt[2]
    d['x_fit'] = x
    d['y_fit'] = y_function(x, *popt)
    d['r_squared'] = r2_score(y, d['y_fit'])
    
    print('popt: %g  %g  %g' % (d['a'], d['b'], d['c']))
    
    return d
            
        
    

def fit_exponential_decay(x, y):
    """ Fits exponential decay with a single exponential of form y = offset + amp*exp(-k*x) to y data """    

    def y_single_exp(x_data, offset, amp, k):
        print('%g   %g   %g' % (offset, amp, k))
        y_data = offset + np.abs(amp) * (1 - np.exp(-np.abs(k) * x_data))
        return y_data

    min_bounds = [-np.inf, 0, 0.0]
    max_bounds = [np.inf, np.inf, np.inf]
    
    st = [y[-1], y[0] - y[-1], 1.0 / (x[-1] - x[0])]
    st = [1000, 100, 10]
    print(st)
    
    min_bounds = [0, 0, 0]
    max_bounds = [10000, 10000, 10000]
    
    try:
        print('hello')
        
        popt, pcov, a, b, c = curve_fit(y_single_exp, x, y, st,
                                        bounds=(min_bounds, max_bounds),
                                        full_output=True)
        
        print('popt')
        print(popt)
        
    except:
        print('fit exponential decay failed - setting decay rate to NaN')
        popt = [y[-1], y[0]-y[-1], np.nan]
        exit(1)

    d = dict()
    d['offset'] = popt[0]
    d['amp'] = popt[1]
    d['k'] = popt[2]
    d['x_fit'] = x
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
    d['x'] = x[:,1]
    d['y'] = y
    d['intercept'] = results.params[0]
    d['slope'] = results.params[1]
    d['x_fit'] = d['x']
    d['y_fit'] = d['intercept'] + d['slope'] * d['x']
    
    return d


def remove_outliers(x, y, poly_deg=2, med_threshold=5):
    """ Tries to remove outliers """
    
    keep_going = True
    
    try:
        while keep_going:
            p_fit = poly.fit(x, y, deg = poly_deg)
            
            y_est = p_fit(x)
            
            abs_y_diff = np.abs(y - y_est)
            med_y_diff = np.median(abs_y_diff)
            
            if (np.max(abs_y_diff) > (med_threshold * med_y_diff)):
                # Points deviate, drop the furthest one
                bi = np.argmax(abs_y_diff)
                
                x = np.delete(x, bi)
                y = np.delete(y, bi)
            else:
                keep_going = False
    except:
        x = x
        y = y
    
    return (x,y)