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

def noise_floor_est(n,scale,trace,o,ch):
    step = int(n/2)
    k = 0
    bins = np.arange(0,len(trace),step)
    xsmooth_seed = np.zeros(len(bins)+1)
    ysmooth_seed = np.zeros(len(bins)+1)

    for i in bins:
    
        nstart = i
        nstop = nstart + n

        yobs = trace[nstart:nstop]
        xobs = np.linspace(nstart,nstop,num=len(yobs))
        mu_data = np.median(yobs)
        sigma_data = np.std(yobs)

    
        if mu_data != 0:
            y = yobs[np.where(yobs<mu_data+(scale*sigma_data))]
            x = xobs[np.where(yobs<mu_data+(scale*sigma_data))]
            p=np.polyfit(x,y,o)

            if i == 0:
                xsmooth_seed[k-1]=min(xobs)
                ysmooth_seed[k-1]=np.polyval(p,min(xobs))
            elif i==max(bins):
                xsmooth_seed[k+1]=max(xobs)
                ysmooth_seed[k+1]=np.polyval(p,max(xobs))
            else:         
                xsmooth_seed[k] = min(xobs) + (max(xobs)-min(xobs))/2
                ysmooth_seed[k] = np.polyval(p,xsmooth_seed[k])
        else:
            xsmooth_seed[k]= min(xobs) + (max(xobs)-min(xobs))/2
            ysmooth_seed[k]= 0.0
        k=k+1
    f = interpolate.interp1d(xsmooth_seed, ysmooth_seed,kind='slinear')
    ynoise = f(ch)
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
    noise,upper= noise_floor_est(300,1,trace_norm,2,ch)
    plt.figure()
    plt.plot(ch,trace_norm,ch,noise,ch,upper)
    gigul.show_peaks(trace,ch,peaks,peak_est,idir+'picks_'+fname)
