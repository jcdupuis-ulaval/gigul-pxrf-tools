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


# File setup for data and results################################
tdir = '../data/Terrain/'
ddir = '../results/CSV/denoised/'
rdir = '../results/peaks/'
vdir = '../results/statistics/'
idir = '../results/PNG/denoised/'
fname_list = os.listdir(path=tdir)
b = len(fname_list)
a = 0
for value in fname_list:
    fname= fname_list[a]
    fname=fname[:-4]
    a = a+1
    # Filter parameters #############################################
    amp_sensitivity=0.1    # Only values above mu + sigma*amp_threshold will be considered for peaks 
    slope_sensitivity=0.05  # Only values below mu - sigma threshold will be considered for peaks
    peak_half_width = 5     # Width of the data selection to fit polynomial function
    #################################################################
    
    flist = os.listdir(path=ddir)
    flist = [i for i in flist if i.startswith(fname)]
    
    for fname in flist:
        data=np.genfromtxt(ddir+fname,delimiter=',')
        # We can try to scale the data between 0 and 1 to make picking easier 
        # with the thresholds
        trace = data[:,1]
        trace_norm = gigul.scale_trace(trace)
        ch = data[:,0] 
        # We generate our peak estimates on the normalized data
        peak_est = gigul.estimate_peaks(trace_norm,ch,amp_sensitivity,slope_sensitivity,fname)
        # We pick on the real data starting from our pick estimates
        peaks = gigul.refine_peaks(trace,peak_est,peak_half_width,ch,rdir+'picks_'+fname)
    
        gigul.show_peaks(trace,ch,peaks,peak_est,idir+'picks_'+fname)
    if a < b:
        continue
    if a == b:
        break
