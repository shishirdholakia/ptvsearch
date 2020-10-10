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
    
    def __init__(self,lc,threshold=0.1,maximum_frequency=100):
        self.lc = lc.remove_nans()
        self.threshold = threshold
        self.maximum_frequency = maximum_frequency
        self.ls = LombScargle(lc.time, lc.flux, normalization='psd')
        self.frequency, self.power = self.ls.autopower()
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

    
    def get_discard_peaks(self,discard_threshold=0.01,discard_all=False,minimum_frequency=10,remove_closer=None):
        
        """
        Gets peaks to discard from the periodogram. By default, it finds peaks
        that are smaller than the main peaks of the delta scuti (the ones that
        you might want to use to fit for PTVs) but still large enough to cause unwanted
        noise in the final estimate of PTVs. Can also be used to find all the peaks
        if discard_all=True, in which case it will yield all the peaks above a threshold.
        This is useful for getting rid of all of the peaks in the lightcurve,
        for instance if you wish to search for transits.
        
        """
        
        if remove_closer==None:
            remove_closer = np.inf
        self.discard_peaks = []
        if discard_all:
            discard_peaks, properties = find_peaks(self.power, prominence=[discard_threshold,1])
        else:
            discard_peaks, properties = find_peaks(self.power, prominence=[discard_threshold,self.threshold])
        for peak in discard_peaks:
            if self.frequency[peak]<self.maximum_frequency and self.frequency[peak]>minimum_frequency:
                a = np.abs(self.frequency[discard_peaks] - self.frequency[peak])
                ma = np.ma.masked_equal(a, 0.0, copy=False)
                if np.any(ma<remove_closer):
                    if self.power[peak] > self.power[discard_peaks[ma.argmin()]]:
                        self.discard_peaks.append(peak)
                        print("Discarded peak at: " + str(self.frequency[discard_peaks[ma.argmin()]]))
                else:    
                    self.discard_peaks.append(peak)
        return self.discard_peaks
    
    def get_discard_peaks_v2(self,discard_threshold=0.01,discard_all=False,minimum_frequency=10,remove_closer=None):
        
        """
        Gets peaks to discard from the periodogram. By default, it finds peaks
        that are smaller than the main peaks of the delta scuti (the ones that
        you might want to use to fit for PTVs) but still large enough to cause unwanted
        noise in the final estimate of PTVs. Can also be used to find all the peaks
        if discard_all=True, in which case it will yield all the peaks above a threshold.
        This is useful for getting rid of all of the peaks in the lightcurve,
        for instance if you wish to search for transits.
        
        """
        
        if remove_closer==None:
            remove_closer = np.inf
        self.discard_peaks = []
        if discard_all:
            discard_peaks, properties = find_peaks(self.power, prominence=[discard_threshold,1])
        else:
            discard_peaks, properties = find_peaks(self.power, prominence=[discard_threshold,self.threshold])
        for peak in discard_peaks:
            if self.frequency[peak]<self.maximum_frequency and self.frequency[peak]>minimum_frequency:
                self.discard_peaks.append(peak)
        #now remove elements if they are within the remove_closer threshold
        if remove_closer is not np.inf:
            discard_peaks = self.discard_peaks
            self.discard_peaks = []
            for peak in discard_peaks:
                highest_peak = True
                for close_peak in discard_peaks:
                    if np.abs(self.frequency[peak]-self.frequency[close_peak]) < remove_closer and self.power[peak] < self.power[close_peak]:
                        highest_peak = False
                        print("Discarded peak at: " + str(self.frequency[close_peak]))
                if highest_peak:
                    self.discard_peaks.append(peak)
        return self.discard_peaks
        
        
    def get_periodogram_quality(self):
        pass
    
    def rank_by_quality(self):
        pass
    
    def get_periodogram_peaks(self):
        return np.array([self.frequency[self.peaks],self.power[self.peaks]])
    
#helper function to create guess parameters from list of peaks    
    def create_guess_values(self,peaks):
        guesses = []
        for index, peak in enumerate(peaks):
            y_fit = self.ls.model(self.lc.time, self.frequency[peak])
            amp = np.sqrt(2)*np.std(y_fit)
            phase = self.lc.time[np.abs(y_fit - 1.0).argmin()]
            guesses.append([self.frequency[peak], amp, phase])
        return guesses
     
    def clean_lc(self, guesses,joint_closer=False, window_size=0.08643*2,amp_tol=None):
        
        """
        Cleans lightcurve by removing sections of the periodogram that aren't
        near the peaks.
        
        Applies fit_sine_simplex one at a time for each peak in the lightcurve
        
        Takes input guess which is array of freqs,ampls,phases for each peak in the periodogram
        
        Also takes parameter joint_closer. 
            If joint_inwindow is True, the algorithm will check if there are
            multiple peaks in a single window. If this is true, it will fit these
            peaks jointly
            
            If joint_closer is None, the algorithm will fit every peak independently.
        
        returns the cleaned lightcurve and an array of the resulting fit parameters
        
        """
        #make a copy of the lc object, subtract the sines from the copy, and return it
        lc = self.lc.copy().remove_nans()
        
        # if there are close peaks, bunch them together into one guess to fit jointly
        #for instance if you have [[1,10,0], [1.1,20,0], [5,25,0]] and window_size=0.5:
        #make that [[1,10,0,1.1,20,0], [5,25,0]]
        if joint_closer:
            guesses_bunched = []
            for freq,amp,phase in guesses:
                highest_peak = True
                guess_bunched = []
                for freq_close,amp_close,phase_close in guesses:
                    if np.abs(freq-freq_close) < window_size/2:
                        if amp < amp_close:
                            highest_peak = False
                        guess_bunched.append(freq_close)
                        guess_bunched.append(amp_close)
                        guess_bunched.append(phase_close)
                if highest_peak:
                    if len(guess_bunched)==3:
                        print("Fitting: " + str(guess_bunched))
                    else:
                        print("Jointly fitting: " + str(guess_bunched))
                    guesses_bunched.append(guess_bunched)
            guesses = guesses_bunched

            
        params_array = []
        for guess_params in guesses:
            params = self.fit_sine_simplex(steps=500, hf_width=window_size/2, window_samples=100, guess=guess_params,amp_tol=amp_tol)
            sinusoid = 0
            for i in range(len(params)//3):
                sinusoid += params[3*i+1]*np.sin(2*params[3*i]*np.pi*(lc.time - params[3*i+2])) #
            lc.flux = lc.flux-sinusoid
            params_array.append(params)
        return lc,params_array
        
    def fit_sine_simplex(self, steps, hf_width, window_samples, guess, lc=None,freq_tol =0.01,amp_tol=None):
        """ 
        Helper function to fit a sine wave to a signal by minimizing the 
        significance of the peak in frequency space (from a Lomb-Scargle)
        
        guess argument should be tuple of size 3n containing:
            frequency initial guess
            amplitude initial guess
            phase initial guess
            frequency initial guess 2 (if joint fitting)
            .
            .
            etc.
            
        freq_tol: amount frequency is allowed to vary in the fit compared to 
        the window size, default 1% of the window size
        """
        if lc is None:
            lc = self.lc
            
        #create a window around the tallest peak
        tallest_peak= max([guess[3*i+1] for i in range(len(guess)//3)])
        tallest_peak_index = guess.index(tallest_peak)-1
            
            
        freq_grid = np.linspace(guess[tallest_peak_index]-hf_width, guess[tallest_peak_index]+hf_width, window_samples)
        
        #create bounds to prevent the freq from straying outside the window
        #and also prevent the phase from exploring a larger space than it needs to
        if amp_tol is None:
            bnds = [
                    [(guess[3*i]-hf_width*freq_tol,guess[3*i]+hf_width*freq_tol),
                     (None,None),
                     (guess[3*i+2]-(1/guess[3*i]),guess[3*i+2]-(1/guess[3*i]))]
                    
                    for i in range(len(guess)//3)
                ]
        else:
            bnds = [
                    [(guess[3*i]-hf_width*freq_tol,guess[3*i]+hf_width*freq_tol),
                     (guess[3*i+1]-amp_tol, guess[3*i+1]+amp_tol),
                     (guess[3*i+2]-(1/guess[3*i]),guess[3*i+2]-(1/guess[3*i]))] #
                    
                
                for i in range(len(guess)//3)
                ]
        bnds = [item for sublist in bnds for item in sublist]
        params = optimize.minimize(fitfunc, guess, args=(lc,freq_grid), method='L-BFGS-B', bounds=bnds)
        return params.x
        
        
        
        
### HELPER FUNCTIONS
def fitfunc(params,*args):
    """
    lc, grid, guess = x
    freq, amp, phase = params
    
    Performs a joint fit to any number of sine waves, where each set of 3
    params correspond to the period, amplitude, and phase of a given sine wave
    """
    assert len(params)%3 == 0
    lc,grid = args
    sinusoid = 0
    for i in range(len(params)//3):
        sinusoid += params[3*i+1]*np.sin(2*params[3*i]*np.pi*(lc.time - params[3*i+2])) #
    ls = LombScargle(lc.time,
                     lc.flux-sinusoid)
    power = ls.power(grid, method="cython")
    return power.max()

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

    
    