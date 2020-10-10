#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 02:16:31 2019

@author: shashank
"""
import os
import pandas as pd
import lightkurve as lk
from del_scu_search import Periodogram
import matplotlib.pyplot as plt
from ptv_search import find_phase_OC, plot_LAT

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
            
            #remove nans from lc, save an LC figure
            lc = lc.remove_nans()
            lc = lc.remove_outliers(sigma_upper=5, sigma_lower=20)
            ax = lc.plot()
            ax.figure.savefig('lightcurves/'+str(ticid)+'_lc.png')
            plt.close()
            
            #create a lob-scargle periodogram to check if it's a delta scuti
            periodogram = Periodogram(lc,threshold=0.01)
            del_scu = periodogram.is_delta_scuti()
            print("Delta Scuti?: " + str(del_scu))
            
            #if it's a delta scuti, save a figure of the periodogram with the peaks highlighted
            if del_scu:
                frequency,power,peaks = periodogram.frequency, periodogram.power,periodogram.peaks
                peak_freqs = frequency[peaks]
                print(peak_freqs)
                freqs = []
                for i in peak_freqs:
                    pg = lc.to_periodogram(minimum_frequency = i-0.3,maximum_frequency = i+0.3, oversample_factor = 500, nyquist_factor = 4)
                    freqs.append(pg.frequency_at_max_power.value)
                print(freqs)
                discard_peaks = periodogram.get_discard_peaks_v2(minimum_frequency=5,discard_threshold=0.001,discard_all=True,remove_closer=0.1)
                
                guess = periodogram.create_guess_values(discard_peaks)
                
                plt.clf()
                lc_clean, params = periodogram.clean_lc(guess)
                ax = lc_clean.plot()
                ax.figure.savefig('clean lightcurves/'+str(ticid)+'_lc.png')
                plt.close()
                df = pd.DataFrame(data=params, columns = ['Frequency','Amplitudes','Phase'])
                df.to_csv('frequencies/'+str(ticid)+'_pg.csv')
                
                
                plt.clf()
                fig, axs = plt.subplots(2,sharex=True, sharey=True)
                fig.suptitle('Periodograms')
                axs[0].plot(frequency,power,c='blue',zorder=3)
                axs[0].scatter(frequency[peaks],power[peaks],c='black',zorder=2)
                axs[0].scatter(frequency[discard_peaks],power[discard_peaks],c='red',zorder=1)
                
                periodogram_clean = Periodogram(lc_clean,maximum_frequency=100, threshold = 0.01)
                frequency_clean,power_clean,peaks_clean = periodogram_clean.frequency, periodogram_clean.power,periodogram_clean.peaks
                axs[1].plot(frequency_clean,power_clean,c='blue',zorder=3)
                axs[1].scatter(frequency_clean[peaks_clean],power_clean[peaks_clean],c='black',zorder=2)
                discard_peaks_clean = periodogram_clean.get_discard_peaks_v2(minimum_frequency=5,discard_threshold=0.001,discard_all=True,remove_closer=0.1)
                axs[1].scatter(frequency_clean[discard_peaks_clean],power_clean[discard_peaks_clean],c='red',zorder=1)                

                
                axs[0].set_xlabel('Frequency (cycles/day)')
                axs[0].set_ylabel('Power')
                axs[0].set_xlim([0,100])
                axs[1].set_xlabel('Frequency (cycles/day)')
                axs[1].set_ylabel('Power')
                axs[1].set_xlim([0,100])
                plt.savefig('periodograms/'+str(ticid)+'_pg.png',dpi=200)
                plt.close()
                
                try:
                    #plot the OC plot for the delta scuti
                    oc_modes = find_phase_OC(lc,freqs)
                    plot_LAT('o-c plots/'+str(ticid)+'_o_c.png', oc_modes)
                except:
                    print("Failed to make LAT plot for TIC " + str(ticid))
        except():
            print("No LCs found. Skipping")

            

            
        plt.close('all')
        