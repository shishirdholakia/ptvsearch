#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 01:03:15 2019

@author: shashank
"""
import lightkurve as lk
import matplotlib.pyplot as plt
import scipy
pi = np.pi    
def find_phase_OC(lc):
    #Download light curve for each sector and stitch it into one light curve
    starid =
    frequency_list = []
    
    #iterate over all the pulsation modes of the star find the OC
    for freq in frequency_list:
        def sinfunc(t,c, A):  return A * np.sin(2*pi*freq*(t+c)) + 1.0
    
    
    
        #subdivide every light curve into 100 sections
        num_sections = 100
        time = list(divide_chunks(lc.time,num_sections))
        flux = list(divide_chunks(lc.flux,num_sections))
        import numpy as np
        periodlist = []
        mediantimelist = []
        #iterate over every window
        for index, interval in enumerate(time):
            #make the window into a TESSLightCurve object
            lc_chunk = lk.LightCurve(interval, flux[index])
            guess = np.array([guess_phase,guess_amp])
            popt, pcov = scipy.optimize.curve_fit(sinfunc, interval, flux[index], p0=guess)
            periodlist.append(popt[0])
            mediantimelist.append(np.median(interval))
        return(periodlist,mediantimelist)
def divide_chunks(l, n): 
  
# looping till length l 
for i in range(0, len(l), n):  
    yield l[i:i + n]
    
def plot_LAT(path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    arr = scipy.stats.binned_statistic(mediantimelist, periodlist, 'mean', bins=20)
    ax.scatter(mediantimelist,periodlist,s=5)
    ax.scatter(arr[1][0:-1],arr[0],c='red',s = 20)
    plt.savefig(path)

def find_period_OC():
    pass

def PTV_model():
    pass

def do_joint_fit():
    pass

