#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 00:51:15 2019

@author: shashank
"""
import numpy as np
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
    
    def get_periodogram_peaks(self):
        return np.array([self.frequency[self.peaks],self.power[self.peaks]])
    
    def clean_lc(self):
        """
        Cleans lightcurve by removing sections of the periodogram that aren't
        near the peaks.
        """
        
        
        return 
def sine(x,t,b,a):
    return a*np.sin(2*(1/t)*np.pi*(x - b))
    
class SinModel:
    
    """
    Base class for a variety of models of a Delta Scuti star's
    pulsations, ranging from a simple sine function to a (not yet implemented)
    full physical model with light arrival time modulations to the phase and the
    frequency
    """
    
    def __init__(self,t,freq, phase, amp,flux=None):
        
        if (type(freq) is float or type(freq) is int) \
        and (type(phase) is float or type(phase) is int) \
        and (type(amp) is float or type(amp) is int):
            
            self.freqs = [freq]
            self.phases = [phase]
            self.amps = [amp]
        else:
            self.freqs = freq
            self.phases = phase
            self.amps = amp          
        
        self.time = np.array(t)
        if flux is None:
            self.flux = np.zeros(len(self.time))
            for f,p,a in zip(self.freqs,self.phases,self.amps):
                self.flux += SinModel.model(self,f, p, a)
        else:
            self.flux = flux
    
    def model(self,freq, phase, amp):

        return amp * np.sin(2*np.pi*freq*(self.time + phase)) + 1.0
    
    def add_noise_model(self,sigma,dist="gaussian"):
        self.flux = self.flux + np.random.normal(0,sigma,len(self.time))
        return self.flux
    
    def __add__(self, other):
        assert len(self.time)==len(other.time)
        assert len(self.flux) == len(other.flux)
        time  = self.time
        total_flux = self.flux + other.flux
        freqs = self.freqs+other.freqs
        phases = self.phases+other.phases
        amps = self.amps + other.amps

        return self.__class__(time,freqs,phases,amps,flux=total_flux)
        
    
class PM_Model(SinModel):
        def __init__(self, t,freq, phase, A, pm_period=0,pm_phase=0,pm_A=0,flux = None):
            super().__init__(t,freq, phase, A)
            print(self.freqs,self.amps,self.phases,self.flux)
            self.pm_period = pm_period
            self.pm_phase = pm_phase
            self.pm_A = pm_A
            if flux is None:
                self.flux = self.model(pm_period,pm_phase,pm_A)
            else:
                self.flux = flux
            
        def model(self,pm_period,pm_phase,pm_A):
            """
            Same as simple sine function, but uses a simple sine as a time
            varying function to describe amplitude
            """
            phase_func = sine(self.time,self.pm_period,self.pm_phase,self.pm_A)
            flux = np.zeros(len(self.time))
            for freq,amp,phase in zip(self.freqs,self.amps,self.phases):
                flux += amp * np.sin(2*np.pi*freq*(self.time + (phase_func+phase))) + 1.0
            return flux
    
class FM_model(SinModel):
    pass

    
    