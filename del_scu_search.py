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
    
    def __init__(lc,threshold=0.1)
        self.lc = lc
        self.threshold
        self.frequency, self.power = LombScargle(lc.time, lc.flux).autopower()
        self.peaks, properties = find_peaks(power, prominence=threshold)
        
    def is_delta_scuti():

        if len(self.peaks)>0:
            return [self.frequency,self.power,self.peaks]
        elif len(self.peaks)==0:
            return False
        
    def get_periodogram_quality():
        pass
    
    def rank_by_quality():
        pass
    
    def get_periodogram_peaks():
        return np.array([self.frequency[self.peaks],self.power[self.peaks]])
    
    def clean_lc():
        """
        Cleans lightcurve by removing sections of the periodogram that aren't
        near the peaks.
        """
        return 