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
#fname ='TDCAEXD332727z_paire5'
ddir = '../results/CSV/denoised/'
rdir = '../results/peaks/'
idir = '../results/PNG/'
# Filter parameters #############################################
amp_threshold=3e-2   # Only values above this thershold will be considered for peaks 
slope_threshold =-2e-4  # Only values below this threshold will be considered for peaks
peak_half_width = 5     # Width of the data selection to fit polynomial function
#################################################################

flist = os.listdir(path=ddir)
    
for fname in flist:
    data=np.genfromtxt(ddir+fname,delimiter=',')
    # We can try to scale the data between 0 and 1 to make picking easier 
    # with the thresholds
    trace_norm = (data[:,1]-data[:,1].min())/(data[:,1].max()-data[:,1].min())
    trace = data[:,1]
    ch = data[:,0] 
    # We generate our peak estimates on the normalized data
    peak_est = gigul.estimate_peaks(trace_norm,amp_threshold,slope_threshold)
    # We pick on the real data starting from our pick estimates
    peaks = gigul.refine_peaks(trace,peak_est,peak_half_width,ch,rdir+'picks_'+fname)
    gigul.show_peaks(trace,ch,peaks,peak_est,idir+'picks_'+fname)
