#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 00:51:15 2019

@author: shashank
"""
from astropy.stats import LombScargle
from scipy.signal import find_peaks
import scipy.fftpack

class Periodogram:
    
    def __init__(self,lc,threshold=0.1):
        self.lc = lc
        self.threshold = threshold
        self.frequency, self.power = LombScargle(lc.time, lc.flux).autopower()
        self.peaks, properties = find_peaks(self.power, prominence=threshold)
        
        
    def is_delta_scuti(self):

        if len(self.peaks)>0:
            return True
        elif len(self.peaks)==0:
            return False

    
    def get_discard_peaks(self,discard_threshold=0.01):
        self.discard_peaks, properties = find_peaks(self.power, prominence=[discard_threshold,self.threshold])
        return self.discard_peaks
        
        
    def get_periodogram_quality(self):
        pass
    
    def rank_by_quality(self):
        pass
    
    def get_periodogram_peaks():
        return np.array([self.frequency[self.peaks],self.power[self.peaks]])
    
    def clean_lc(self):
        """
        Cleans lightcurve by removing sections of the periodogram that aren't
        near the peaks.
        """
        
        
        return 
    
class SinModel:
    
    """
    Base class for a variety of models of a Delta Scuti star's
    pulsations, ranging from a simple sine function to a (not yet implemented)
    full physical model with light arrival time modulations to the phase and the
    frequency
    """
    
    def __init__(self,t,freq, phase, A):
        self.freq = freq
        self.phase = phase
        self.A = A
        
        self.time = t
        self.flux = self.model(t,freq, phase, A)
    
    def model(self):

        return self.A * np.sin(2*np.pi*self.freq*(self.time + self.phase)) + 1.0
    
    def noise_model(sigma,type="gaussian"):
        if type=="gaussian":
            return self.model + np.random.normal(0,sigma,len(self.time))
        
    
class PM_Model(SinModel):
    
        def model(t,freq, c, A, pm_freq,pm_c,pm_A):
        """Same as simple sine function, but uses a simple sine as a time
        varying function to describe amplitude"""
        return A * np.sin(2*pi*freq*(t + self.simple_sin(t,pm_freq,pm_c,pm_A))) + 1.0
class 

    
    