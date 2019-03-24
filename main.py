#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 02:16:31 2019

@author: shashank
"""
import os
import pandas as pd
import lightkurve as lk
from del_scu_search import is_delta_scuti
import matplotlib.pyplot as plt
from ptv_search import *

if __name__=='__main__':
    

    if not os.path.exists('lightcurves/'):
        os.makedirs('lightcurves/')
    if not os.path.exists('periodograms/'):
        os.makedirs('periodograms/')
    if not os.path.exists('o-c plots/'):
        os.makedirs('o-c plots/')

    astars = pd.read_csv('tessAstars.csv')
    for star in astars['ID']:
        ticid = 'TIC ' + str(star)
        print(ticid)
        try:
            lcs = lk.search_lightcurvefile(ticid).download_all()
            lc = None
            
            #concatenate the available TESS sectors crudely
            for lcfile in lcs:
                if lc is None:
                    lc = lcfile.PDCSAP_FLUX.normalize().flatten(window_length=201,break_tolerance=10)
                else:
                    lc = lc.append(lcfile.PDCSAP_FLUX.normalize().flatten(window_length=201,break_tolerance=10))
                    
            lc = lc.remove_nans()
            ax = lc.plot()
            ax.figure.savefig('lightcurves/'+str(ticid)+'_lc.png')
            plt.close()
            
            del_scu = is_delta_scuti(lc)
            if del_scu is not False:
                print(del_scu)
                frequency,power,peaks = del_scu
                peak_freqs = frequency[peaks]
                pg = lc.to_periodogram(min_frequency = peak_freqs[0]-0.3,max_frequency = peak_freqs[0]+0.3, oversample_factor = 500, nyquist_factor = 4)
                freqs = [pg.frequency_at_max_power.value]
                print(pg.frequency_at_max_power.value)
                
                plt.clf()
                plt.plot(frequency,power,c='blue')
                plt.scatter(frequency[peaks],power[peaks],c='black')
                plt.xlabel('Frequency (cycles/day)')
                plt.ylabel('Power')
                plt.xlim([0,40])
                plt.savefig('periodograms/'+str(ticid)+'_pg.png')
                plt.close()
                
                periodlist, mediantimelist = find_phase_OC(lc,freqs)
                plot_LAT('o-c plots/'+str(ticid)+'_o_c.png', mediantimelist, periodlist)
        except(TypeError):
            print("No LCs found" Skipping)

            

            
        plt.close('all')
        