#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 00:51:15 2019

@author: shashank
"""
from astropy.stats import LombScargle
from scipy.signal import find_peaks
import scipy.fftpack

def is_delta_scuti(lc):
    threshold=0.1 #need to experiment more
    frequency, power = LombScargle(lc.time, lc.flux).autopower()
    peaks, properties = find_peaks(power, prominence=threshold)
    if len(peaks)>0:
        return [frequency,power,peaks]
    elif len(peaks)==0:
        return False
    
def get_periodogram_quality():
    pass

def rank_by_quality():
    pass

def clean_lc(lc):
    """
    Cleans lightcurve by removing sections of the periodogram that aren't
    near the peaks.
    """
    return lc