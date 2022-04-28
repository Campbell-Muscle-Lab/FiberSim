# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 15:22:06 2022

@author: sako231
"""



import os 
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': "Arial"})
plt.rcParams.update({'font.size': 14}) 

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

    
ROOT = os.path.dirname(__file__)

### fit functions

def fit_straight_line(x, y):
    """ Fits a straight line to data """

    import statsmodels.api as sm
    
    x_data = x
    y_data = y

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

    
    x_fit = np.linspace(x_data[0], x_data[-1], 1000)
    y_fit = d['slope'] * x_fit + d['intercept'] 
    
    
    d['x_fit'] = x_fit
    d['y_fit'] = y_fit

    return d

def fit_exponential_recovery(x, y, n=1):
    """ Fits exponential recovery with a single exponential of form y = offset + amp*(1 - exp(-k*x)) to y data """
    
    def y_single_exp(x_data, offset, amp, k):
        y = np.zeros(len(x_data))
        for i,x in enumerate(x_data):
            y[i] = offset + amp*np.exp(-k*(x))
        return y   

    min_bounds = [0.0, - np.inf, 0.0]
    max_bounds = [np.inf, np.inf, 10.0]
         
    try:
        
        popt, pcov = curve_fit(y_single_exp, x, y, [y[-1], y[0]-y[-1], -np.log(0.5)/0.5*(np.amax(x)+np.amin(x))],
                               maxfev=10000,
                               bounds=(min_bounds, max_bounds))
        
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

### Prepare figure

fig = plt.figure(constrained_layout=True, figsize=(6, 6), dpi = 200)
spec = fig.add_gridspec(2, 1, height_ratios = [1,1])

ax = []

ax0 = fig.add_subplot(spec[0,0]) # Force traces
ax1 = fig.add_subplot(spec[1,0]) # Length traces

ax = [ax0, ax1]

    
for j in range(2):
    
    ax[j].spines['top'].set_visible(False)
    ax[j].spines['right'].set_visible(False)
    
    for axis in ['top','bottom','left','right']:
        ax[j].spines[axis].set_linewidth(1.5)
        
    for tick in ax[j].yaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
        
    for tick in ax[j].xaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
        
    ax[j].tick_params(direction = "out", length = 6, width = 1.5)
    ax[j].set_xticks([0,0.25])
    ax[j].spines['bottom'].set_bounds(0,0.25)
    
    ax[j].set_yticks([800, 1100])
    ax[j].spines['left'].set_bounds(800,1100)
    
    ax[j].set_ylabel("HS length\n(nm)", rotation = 0)

ax[1].spines['bottom'].set_visible(True)  

### Prepare exp function to plot

x = np.arange(0,0.25,0.001)

y_exp = 300 + 800*np.exp(-2*x)

ax[0].plot(x,y_exp) 
ax[0].set_title("'length_fit_mode' set to 'exponential' (default)")


ax[1].plot(x,y_exp) 
ax[1].set_title("'length_fit_mode' set to 'linear'")
ax[1].set_xlabel("Time (s)")

# ax[1].set_ylim([850,1150])            
# ax[1].set_xlim(1.19,1.30)

# color_set = [ "#a0001bff", "#00297ffd", "#7ea8beff", "#df2935b3"]

# label_leg = ["model_1", "model_2", "model_3"]
    
    # ### Fitting
    
    # time_interval = [0.45, 0.8]
    
    # time_fit = [0.50, 0.8]
    

    # for file in os.listdir(data_folder):
    #     if file.endswith('.txt'):
        
    #         data_file_string = os.path.join(data_folder, file)
    
    #         pd_data = pd.read_csv(data_file_string, delimiter='\t')       
        
    #         d_fit = pd_data.loc[(pd_data['time'] >= time_fit[0]) &
    #           (pd_data['time'] <= time_fit[-1])]
        
    #         time_data = d_fit["time"] - time_fit[0]
            
    #         # Plot force
            
    #         ax[0].plot('time', 'force', data = pd_data, color = color_set[i], linestyle = "-", label = label_leg[i]) 
    #         ax[0].legend()
            
    #         # Plot length
        
    #         # Fitting
            
    #         fitted_curve = fit_straight_line(time_data.to_numpy(), d_fit["hs_length"].to_numpy())
            
    #         fitted_curve_exp = fit_exponential_recovery(time_data.to_numpy(), d_fit["hs_length"].to_numpy())
            
    #         vmax = -fitted_curve['slope']/pd_data["hs_length"].loc[0] # velcoity in ML/s
            
    #         vmax_exp = fitted_curve_exp['amp'] * fitted_curve_exp['k']/pd_data["hs_length"].loc[0]
            
    #         ax[1].plot('time', 'hs_length', data = pd_data, color = color_set[i], linestyle = "-", label =  "V_max = %0.2f ML/s" % vmax_exp) 
            
    #         ax[1].plot(fitted_curve_exp['x_fit'] + time_fit[0], fitted_curve_exp['y_fit'], color = 'k', linewidth = 1.5, linestyle = "--") 
            
    #         ax[1].set_ylabel("HS length (nm)", rotation = 90, labelpad = 20)
            
    #         ax[1].set_xlabel("Time (s)", labelpad = 10)
            
    #         ax[1].legend()
    #         #ax[1].set_ylim([850,1150])            
    #         #ax[1].set_xlim(1.19,1.30)

            

    
    
    

    
    
    
    



