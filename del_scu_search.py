#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 00:51:15 2019

@author: shashank
"""
import numpy as np
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import scipy.fftpack
from scipy import optimize

class Periodogram:
    
    def __init__(self,lc,threshold=0.1,maximum_frequency=40):
        self.lc = lc
        self.threshold = threshold
        self.maximum_frequency = maximum_frequency
        self.frequency, self.power = LombScargle(lc.time, lc.flux,normalization='psd').autopower()
        self.discard_peaks = []
        peaks, properties = find_peaks(self.power, prominence=threshold,distance=500)
        
        self.peaks = []
        for peak in peaks:
            if self.frequency[peak]<self.maximum_frequency:
                self.peaks.append(peak)
        
    def is_delta_scuti(self):

        if len(self.peaks)>0:
            return True
        elif len(self.peaks)==0:
            return False

    
    def get_discard_peaks(self,discard_threshold=0.01):
        discard_peaks, properties = find_peaks(self.power, prominence=[discard_threshold,self.threshold])
        for peak in discard_peaks:
            if self.frequency[peak]<self.maximum_frequency:
                self.discard_peaks.append(peak)        
        return self.discard_peaks
        
        
    def get_periodogram_quality(self):
        pass
    
    def rank_by_quality(self):
        pass
    
    def get_periodogram_peaks(self):
        return np.array([self.frequency[self.peaks],self.power[self.peaks]])
    
    def clean_lc(self, peaks=self.peaks, ampls, phases):
        
        """
        Cleans lightcurve by removing sections of the periodogram that aren't
        near the peaks.
        
        Applies remove_sine_simplex one at a time for each peak in the lightcurve
        
        """
        #for i, peak in enumerate(peaks):
            #fit_sine_simplex()
        pass
        
        
    def fit_sine_simplex(self, lc=self.lc, freq_sf, amp_sc, phase_sc, steps, hf_width, window_samples, guess):
        """ 
        Helper function to fit a sine wave to a signal by minimizing the 
        significance of the peak in frequency space (from a Lomb-Scargle)
        
        guess argument should be tuple containing:
            frequency initial guess
            amplitude initial guess
            phase initial guess
        """
        #calculate initial periodogram significance
        freq_grid = np.linspace(guess[0]-hf_width, guess[0]+hf_width, window_samples)
        power = LombScargle(lc.time, lc.flux, lc.flux_err).power(freq_grid, method="cython")
        initial_sig = 1-ls.false_alarm_probability(power.max())
        
        params = minimize(fitfunc, guess, args=(lc,freq_grid), method='Nelder-Mead')
        
        return params[1]*np.sin(2*params[0]*np.pi*(lc.time - params[2]))
        
        
        
        

def fitfunc(params,*args):
    """
    lc, grid, guess = x
    freq, amp, phase = params
    """
    
    lc,grid = args
    ls = LombScargle(lc.time,
                     lc.flux-params[1]*np.sin(2*params[0]*np.pi*(lc.time - params[2])),
                     lc.flux_err)
    power = ls.power(grid, method="cython")
    return 1-ls.false_alarm_probability(power.max())

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
            varying function to describe phase
            """
            phase_func = sine(self.time,pm_period,pm_phase,pm_A)
            flux = np.zeros(len(self.time))+1.0
            for freq,amp,phase in zip(self.freqs,self.amps,self.phases):
                flux += amp * np.sin(2*np.pi*freq*(self.time - (phase_func+phase)))
            return flux
    
class FM_model(SinModel):
    pass

    
    