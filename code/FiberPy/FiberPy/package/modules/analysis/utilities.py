# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:20:25 2020

@author: kscamp3
"""

import numpy as np
import os

def fit_pCa_data(x,y):
    """ Fits Hill-curve to x-y data """
    
    from scipy.optimize import curve_fit

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

def save_figure_to_file(f, im_file_string, dpi=250, verbose=1):
    """ Writes an image to file """
    
    dir_path = os.path.dirname(im_file_string)
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
        
    if (verbose):
        print('Saving figure to: %s' % im_file_string)
        
    f.savefig(im_file_string, dpi=dpi)