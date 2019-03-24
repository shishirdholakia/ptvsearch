#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 02:16:31 2019

@author: shashank
"""
import pandas as pd
import lightkurve as lk
from del_scu_search import is_delta_scuti
import matplotlib.pyplot as plt
from ptv_search import *

if __name__=='__main__':
    astars = pd.read_csv('tessAstars.csv')
    for star in astars['ID'][0:1]:
        ticid = 'TIC ' + str(star)
        print(ticid)
        lcs = lk.search_lightcurvefile(ticid).download_all()
        lc = None
        
        #concatenate the available TESS sectors crudely
        for lcfile in lcs:
            if lc is None:
                lc = lcfile.PDCSAP_FLUX
            else:
                lc.append(lcfile.PDCSAP_FLUX)
                
        lc = lc.remove_nans()
        lc.flatten(window_length=101,break_tolerance=50)
        
        del_scu = is_delta_scuti(lc)
        if del_scu is not False:
            pg = lc.to_periodogram(min_frequency = del_scu[0]-0.3,max_frequency = del_scu[0]+0.3, oversample_factor = 500, nyquist_factor = 4)
            freqs = [pg.frequency_at_max_power.value]
            
            periodlist, mediantimelist = find_phase_OC(lc,freqs)
            plot_LAT('o_c.png', mediantimelist, periodlist)
        
        
    plt.show()