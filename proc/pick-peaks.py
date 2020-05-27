# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:35:36 2020
This function looks into the denoised folder for spectra that no longer
have a contaminating background noise. It makes a first estimate of the peaks 
using an amplitud threshold and a slope thershold. The data is then subsampled
and a fourth order polynomial function is fitted to obtain a more accurate 
pick. 
@author: chdup58
"""

import numpy as np
import gigul_pxrf_tools as gigul
import os
from scipy import interpolate
import matplotlib.pyplot as plt

# File setup for data and results################################
#fname ='TDCAEXD332727z_paire5'
ddir = '../results/CSV/denoised/'
rdir = '../results/peaks/'
idir = '../results/PNG/'
# Filter parameters #############################################
amp_sensitivity=0.05    # Only values above mu + sigma*amp_threshold will be considered for peaks 
slope_sensitivity=0.05  # Only values below mu - sigma threshold will be considered for peaks
peak_half_width = 5     # Width of the data selection to fit polynomial function
#################################################################

flist = os.listdir(path=ddir)

def noise_floor_est(n,trace,ch,ch_start):
    ''' To make a noise estimate, we need to determine what is signal and what is noise. We can do this by 
    following a similar approach to what we've done for the background estimates'''
    n=100
    background = gigul.estimate_background(n,scale,trace,o,ch)
    return ynoise, ynoise+(75.0*scale*sigma_data)


    
for fname in flist:
    data=np.genfromtxt(ddir+fname,delimiter=',')
    # We can try to scale the data between 0 and 1 to make picking easier 
    # with the thresholds
    trace = data[:,1]
    trace_norm = (data[:,1]-(abs(data[:,1]).min()))/(abs(data[:,1]).max()-abs(data[:,1]).min())
    ch = data[:,0] 
    # We generate our peak estimates on the normalized data
    peak_est = gigul.estimate_peaks(trace_norm,amp_sensitivity,slope_sensitivity)
    # We pick on the real data starting from our pick estimates
    peaks = gigul.refine_peaks(trace,peak_est,peak_half_width,ch,rdir+'picks_'+fname)

    gigul.show_peaks(trace,ch,peaks,peak_est,idir+'picks_'+fname)
