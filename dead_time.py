#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 10:19:05 2022

@author: beriksso
"""

import sys
sys.path.insert(0, '/home/beriksso/TOFu/analysis/benjamin/github/TOFu/functions')
import tofu_functions as dfs
import useful_defs as udfs
import matplotlib.pyplot as plt
import numpy as np


def get_data(shot_number, detector):
    """Return time stamps for tofor/tofu for given shot number and detector."""
    tofu_t = dfs.get_times(shot_number, detector_name=detector)
    tofor_t = udfs.get_tofor_times(shot_number)
    
    return tofu_t, tofor_t[detector]


def plot_delta_t(tofu_t, tofor_t, shot_number):
    """Plot the time difference between time stamps for tofu and tofor."""
    plt.figure('delta t')
    
    # Calculate delta t
    tofu_dt = np.diff(tofu_t)
    tofor_dt = np.diff(tofor_t)

    # Create histogram of delta t  
    width = 2
    bin_edges = np.arange(0, 260+width, width)
    bin_centres = bin_edges[1:] - np.diff(bin_edges)/2
    tofu_h, _ = np.histogram(tofu_dt, bin_edges)
    tofor_h, _ = np.histogram(tofor_dt, bin_edges)
    
    # Plot and configure 
    plt.step(bin_centres, tofu_h, color='k', linestyle='-', label='TOFu')
    plt.step(bin_centres, tofor_h, color='C0', linestyle='--', 
             label='Original DAQ')
    plt.xlabel('$\Delta t_{S1}$ (ns)')
    plt.ylabel('counts')
    plt.xlim(0, 250)
    plt.legend()
    plt.title(f'JPN {shot_number}', loc='left')
    
shot_number = 95776
detectors = dfs.get_dictionaries('S1')
#detectors = dfs.get_dictionaries('S2')
tofu_t = np.array([])
tofor_t = np.array([])
for detector in detectors:
    data = get_data(shot_number, detector)
    tofu_t = np.append(tofu_t, data[0])
    tofor_t = np.append(tofor_t, data[1])
plot_delta_t(tofu_t, tofor_t, shot_number)
