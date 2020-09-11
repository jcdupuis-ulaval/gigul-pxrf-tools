# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:44:00 2020
This set of tools are used to process and pXRF data. 
@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def scale_trace(trace):
    scaled_trace = (trace-(abs(trace).min()))/(abs(trace).max()-abs(trace).min())
    return scaled_trace


def estimate_background(n,scale,trace,o,ch):
    '''
The function remove_background estimates the Bremsstrahlung radiation in the spectra
and removes it from the outpu trace.
# Filter parameters #############################################
ns                  # Width of the window for noise estimate
scale               # SNR Threshold
o                   # Order of the noise approximation 
ch                  # The channels associated with each spectra reading
#################################################################
'''
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

    return ynoise

def remove_background(n,scale,trace,o,fname,ch):
    '''
The function remove_background estimates the Bremsstrahlung radiation in the spectra
and removes it from the outpu trace.
# Filter parameters #############################################
ns                  # Width of the window for noise estimate
scale               # SNR Threshold
o                   # Order of the noise approximation 
fname               # The filename where the denoised data will be stored
ch                  # The channels associated with each spectra reading
#################################################################
'''
    background_estimate = estimate_background(n,scale,trace,o,ch)
    print('Saving denoised trace to file : ' + fname+'.csv')
    np.savetxt(fname+'.csv',np.transpose([ch,trace-background_estimate]),delimiter=',')
    return background_estimate, trace-background_estimate

def calc_amp_threshold(trace,sigma):
    mu_trace=np.median(trace)
    std_trace = np.std(trace)*sigma
    return mu_trace+std_trace

def calc_slope_threshold(trace,sigma):
    mu_trace=np.median(trace)
    std_trace = np.std(trace)*sigma
    return mu_trace-std_trace


def show_clean_trace (ch,trace,ynoise,trace_clean,fname):
    plt.figure()
    plt.plot(3,1,1)
    plt.semilogy(ch,trace)
    plt.grid()
    plt.subplot(3,1,2)
    plt.plot(ch,trace,ch,ynoise,'k')
    plt.grid()
    plt.subplot(3,1,3)
    plt.plot(ch,trace_clean)
    plt.grid()
    print ('Saving figure : '+'denoised_'+fname+'.png')
    plt.savefig(fname+'.png', format='png')

def show_peaks(data,ch,peaks,peak_est,fname):
    plt.figure()
    plt.plot(peaks[:,0],peaks[:,1],'*',peak_est[:,0],peak_est[:,1],'+')
    plt.title(fname)
    plt.legend(['final-pick','first-estimate'])
    plt.plot(ch,data,ch,smooth(data))
    plt.xlabel('Channels')
    plt.ylabel('CPS')
    plt.grid()
    plt.show()
    print ('Saving figure : '+fname+'.png')
    plt.savefig(fname+'.png', format='png')

def smooth (data):
    n = len(data)
    sdata=np.zeros(n)
    sdata[0]=data[0]
    sdata[n-1]=data[n-1]
    for i in np.arange(1,n-1):
        sdata[i]=(data[i-1]+(2*data[i])+data[i+1])/4.0
    return sdata

def estimate_peaks (data,ch,amp_sensitivity,slope_sensitivity):
    # Need to change to include the channels 
    #ch = np.linspace(1,len(data),num=len(data)) # Assign channel numbers 
    data_smooth =smooth(smooth(smooth(data)))
    d = np.gradient(data_smooth)
    dd = smooth(smooth(smooth(np.gradient(d)))) # Generate a second derivative to evaluate slope at zero crossing
   
    #Plot to troubleshoot settings
    #The amplitud threshold and slope threshold values are shown by red lines on the plot
    amp_threshold = calc_amp_threshold(data,amp_sensitivity)
    slope_threshold = calc_slope_threshold(dd,slope_sensitivity)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(ch, data, 'g-')
    ax1.plot(ch,np.ones(len(ch))*amp_threshold,'-r')
    ax2.plot(ch, dd, 'b-')
    ax2.plot(ch,np.ones(len(ch))*slope_threshold,'-r')
    ax1.set_xlabel('ch')
    ax1.set_ylabel('data', color='g')
    ax2.set_ylabel('1st deriv', color='b')
    plt.grid()
    plt.show()
    
    peak=np.zeros((len(d),2))

    for i in np.arange(0,len(d)-1):
        if np.sign(d[i])>np.sign(d[i+1]):
            peak[i,0]=ch[i]
            peak[i,1]=data_smooth[i]
    peak_est = peak[np.where((peak[:,1]>amp_threshold)&(dd<slope_threshold))]
    
    return peak_est

def refine_peaks (data,peak_est,peak_half_width,ch,fname):
    # In order to be able to pick data slices in the spectra, we introduce the concept of an offset
    offset = ch.min()
    ch = ch - offset
    peaks = np.zeros((len(peak_est),2))
    for i in np.arange(0,len(peak_est)):
        pk = int(ch[np.where(ch == (peak_est[i,0]-offset))])
        xstart = pk-peak_half_width
        xstop =  pk+peak_half_width
        if xstart > 0 : # this means that we did not place the start of our data subset on a peak
            x = ch[xstart:xstop]
            y = data[xstart:xstop]

            xnew = np.linspace(xstart,xstop,num=1000)
            p=np.polyfit(x,y,4)
            ynew = np.polyval(p,xnew)
            d = np.gradient(ynew)
            local_peak = np.zeros((1,2))
            for j in np.arange(0,len(d)-1):
                if np.sign(d[j]) > np.sign(d[j+1]):
                    local_peak[0,0]=xnew[j]
                    local_peak[0,1]=ynew[j]
            peaks[i,0]= local_peak[0,0] + offset
            peaks[i,1]= local_peak[0,1] 
            print('Saving refined-picks to file : ' + fname)
            np.savetxt(fname,peaks,delimiter=',')
        elif xstart<0 :
            print ('You cannot sub-divide your dataset on a peak')
            break
    return peaks

# Routine to compute the statistics on the collated peaks  
def calc_stats (raw):
    data=raw[raw[:,1].argsort()] # start by sorting the data to get similar channels next to each other
    m,n = data.shape
    d = np.zeros((m,1))   # Initialize the distance matrix 
    for i in np.arange(m-1): # Compute the distance between adjacent channel numbers needed for the grouping 
        d[i] = data[i+1,1]-data[i,1]

    ch_edges = np.where(d>1)[0] # Define the edges of the data set that should be averaged
    m = len(ch_edges) # Setup the counting variable for going through the list of edges
    # Setup all of the vectors we are going to fill
    mu_ch = np.zeros((m,1))  # Average channel number
    mu_amp = np.zeros((m,1)) # Average amplitude at a given channel
    std_amp = np.zeros((m,1)) # Standard deviation on the amplitude at a given channel
    std_ch = np.zeros((m,1))  # Standard deviation on the position of a given channel
    n = np.zeros((m,1))       # Number of observations that have contributed to the solution
    for i in np.arange(m-1):  # Traverse the list to identify the data subsets to group for analysis
       if i == 0: # It is the first time in this loop 
            ch = data[0:ch_edges[i]+1,1]
            amp = data[0:ch_edges[i]+1,2]
       elif (i>0):
            ch = data[ch_edges[i]+1:ch_edges[i+1]+1,1]
            amp = data[ch_edges[i]+1:ch_edges[i+1]+1,2]
       n[i]=len(ch)          # Note the number of items that contribute to the computed value 
        
       if len(ch)>=1:        # If we have more than one item in our list we can compute the mean and the standard deviation
            mu_ch[i] = np.mean(ch)
            std_ch[i] = np.std(ch)
            mu_amp[i] = np.mean(amp)
            std_amp[i] = np.std(amp)
       elif len(ch)<1:       # If we only have one value we can assign the single value to the list for completeness 
            mu_ch[i] = data[ch_edges[i],1]
            mu_amp[i] =data[ch_edges[i],2]
            n[i]=1
    return mu_ch, std_ch, mu_amp,std_amp,n

################################################################# 
