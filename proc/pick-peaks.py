# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:35:36 2020

@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
import gigul_pxrf_tools as gigul
import os

# File setup for data and results################################
#fname ='TDCAEXD332727z_paire5'
ddir = '../results/CSV/'
rdir = '../results/peaks/'
idir = '../results/PNG/'
# Filter parameters #############################################
flist = os.listdir(path=ddir)
    

for fname in flist:
    
    data=np.genfromtxt(ddir+fname,delimiter=',',skip_header=1)
    ch = np.linspace(1,len(data),num=len(data)) # Assign channel numbers 
    amp_threshold=125.0
    slope_threshold =-10.0
    peak_half_width = 5


    peak_est = gigul.estimate_peaks(data,amp_threshold,slope_threshold)
    peaks = gigul.refine_peaks(data,peak_est,peak_half_width,ch,rdir+'picks_'+fname)
    gigul.show_peaks(data,ch,peaks,peak_est,idir+'picks_'+fname)
