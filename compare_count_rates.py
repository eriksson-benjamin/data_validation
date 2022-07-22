#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 07:53:16 2022

@author: beriksso
"""

'''
Compare the count rate of S1/S2 detectors to fission chambers and TOFOR.
'''

import sys
sys.path.insert(0, '/home/beriksso/TOFu/analysis/benjamin/github/TOFu/functions')
import tofu_functions as dfs
import numpy as np
import useful_defs as udfs
import ppf
import matplotlib.pyplot as plt
udfs.set_nes_plot_style()
def import_tofu(detector, shot_number):
    '''
    Return TOFu time stamps
    '''
    
    # Import time stamps (ns)
    time_stamps = dfs.get_times(shot_number, detector_name=detector)
    
    # Import time offset (ns)
    board, _ = dfs.get_board_name(detector)
    try: offset = dfs.get_offset(board, shot_number)
    except: offset = 0
    
    return time_stamps-offset

def import_tofor(shot_number):
    '''
    Return TOFOR time stamps
    '''
    
    time_stamps = udfs.get_tofor_times(shot_number)
    return time_stamps

def import_kn1(shot_number):
    '''
    Return fission chamber neutron rates
    '''
    
    kn1 = ppf.ppfget(dda='TIN', dtyp='RNT', pulse=shot_number)
    time_centres = kn1[4]
    kn1_rate = kn1[2]
    
    return time_centres, kn1_rate
    

def count_rate(time_stamps, width):
    '''
    Calculate neutron count rate averaged over given time width.
    '''
    
    # x-axis
    time_steps = np.arange(40, 80, width)
    
    # Calculate count rate for given time width
    args = np.searchsorted(time_stamps*1E-9, time_steps)
    count_rate = np.diff(args)/width
    
    # Bin centres
    time_centres = time_steps[1:]-np.diff(time_steps)/2

    return time_centres, count_rate

def plot_count_rate(detector, shot_number, tofor_times, width):
    # Import TOFu time stamps
    tofu_times = import_tofu(detector, shot_number)

    # Calculate count rates
    kn1_bins, kn1_rate = import_kn1(shot_number)
    tofu_bins, tofu_rate = count_rate(tofu_times, width)
    tofor_bins, tofor_rate = count_rate(tofor_times[detector], width)
    
    # Plot count rates
    # ----------------
    plt.figure(detector)
    plt.plot(tofu_bins-40, tofu_rate/tofu_rate.max(), 'k-', label='TOFu')
    plt.plot(tofor_bins-40, tofor_rate/tofor_rate.max(), 'C0-', label='Old DAQ')
    
    # Configure plot
    plt.xlabel('$t_{JET}$ (s)')
    plt.ylabel('$R_n$ $(s^{-1})$')
    plt.xlim(6, 18)
    plt.ylim(0, 1.05)
    plt.title(f'JPN {shot_number}', loc='left')
    plt.title(f'{detector.replace("_", "-")}', loc='right')
    plt.legend()

def main():
    shot_number = 98044
    
    # Import time stamps for all TOFOR detectors
    tofor_times = import_tofor(shot_number)

    plot_count_rate('S1_01', shot_number, tofor_times, 0.006)
    plot_count_rate('S2_01', shot_number, tofor_times, 0.1)
    
if __name__=='__main__':
    main()
    
