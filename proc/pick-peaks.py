# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:35:36 2020

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
amp_threshold=125.0
slope_threshold =-20.0
peak_half_width = 5
#################################################################

flist = os.listdir(path=ddir)
    
for fname in flist:
    data=np.genfromtxt(ddir+fname,delimiter=',')
    trace = data[:,1]
    ch = data[:,0] 
    peak_est = gigul.estimate_peaks(trace,amp_threshold,slope_threshold)
    peaks = gigul.refine_peaks(trace,peak_est,peak_half_width,ch,rdir+'picks_'+fname)
    gigul.show_peaks(trace,ch,peaks,peak_est,idir+'picks_'+fname)
