#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 01:03:15 2019

@author: shashank
"""
import lightkurve as lk
import matplotlib.pyplot as plt
import scipy
import numpy as np
from lmfit import Model 

pi = np.pi  
  
def find_phase_OC(lc,frequency_peaks,guess_phase=0):    
    #iterate over all the pulsation modes of the star find the OC
    for freq in frequency_peaks:
        def sinfunc(t,c, A):  return A * np.sin(2*pi*freq*(t+c)) + 1.0
        #subdivide every light curve into 100 sections
        num_sections = 200
        time = list(divide_chunks(lc.time,num_sections))
        flux = list(divide_chunks(lc.flux,num_sections))
        guess_phase = lc.time[(np.abs(lc.flux[0:1000]-1.0)).argmin()]


        guess_amp = np.std(lc.flux) * 2.**0.5
        bincenter_phases = []
        bincenter_times = []
        #iterate over every window
        for index, interval in enumerate(time):
            #make the window into a TESSLightCurve object
            guess = np.array([guess_phase,guess_amp])
            popt, pcov = scipy.optimize.curve_fit(sinfunc, interval, flux[index], p0=guess, bounds =([guess[0]-0.25/freq,-0.05],[guess[0]+0.25/freq,0.05]))
            bincenter_phases.append(popt[0])
            bincenter_times.append(np.median(interval))
    #fix this to return an array for each frequency
    return(86400*(np.array(bincenter_phases)-np.mean(bincenter_phases)),np.array(bincenter_times))

def sine(x,t,b,a):
    return a*np.sin(2*(1/t)*pi*(x - b))
def line(x,m,c):
    return m*x + c

def ptv_test(bincenter_times,bincenter_phases, t_guess=840,b_guess=450,a_guess=5):
    sinmodel = Model(sine)
    params = sinmodel.make_params(t=t_guess, b=b_guess, a=a_guess)
    y = 86400*(np.array(bincenter_phases)-np.mean(bincenter_phases))
    x = np.array(bincenter_times)
    sin_result = sinmodel.fit(y, params, x=x)
    
    linemodel = Model(line)
    params1 = linemodel.make_params(m=0, c=0)
    line_result = linemodel.fit(y, params1, x=x)
    
    if sin_result.bic-line_result.bic < 0.0:
        return True
    else:
        return False
    
      
def divide_chunks(l, n): 
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n]
    
def plot_LAT(path, bincenter_times, bincenter_phases):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    arr = scipy.stats.binned_statistic(bincenter_times, bincenter_phases, 'mean', bins=5)
    arr1 = scipy.stats.binned_statistic(bincenter_times, bincenter_phases, np.std, bins=5)

    ax.scatter(bincenter_times,bincenter_phases,s=5)
    ax.errorbar(arr[1][0:-1],(np.array(arr[0])-np.mean(bincenter_phases)),yerr = (arr1[0])/np.sqrt(arr[1][1]-arr[1][0]),fmt = 'o', c='red')
    plt.savefig(path)

def find_period_OC():
    pass


def do_joint_fit():
    pass

