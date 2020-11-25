# -*- coding: utf-8 -*-
"""
Created on Mon May 25 15:58:36 2020

@author: debth
"""

import numpy as np
import matplotlib.pyplot as plt
import gigul_pxrf_tools as gigul
import os


# File setup for data and results################################

ddir = '../data/Terrain/'
rdir = '../results/'
fname_list = os.listdir(path=ddir)
b = len(fname_list)
a = 0
for value in fname_list:
    fname= fname_list[a]
    fname=fname[:-4]
    a = a+1
    # Filter parameters #############################################
    ns=50               # Width of the window for noise estimate
    scale = 0.05        # SNR Threshold
    o = 1               # Order of the noise approximation 
    #################################################################
    # load the data file 
    print ('Processing file : '+ddir+fname+'.csv')
    data=np.genfromtxt(ddir+fname+'.csv',delimiter=';',skip_header=1)
    m,n = data.shape
    traces = data[:,1:n]
    nsample,ntraces = traces.shape
    
    
    # Merge data
    if np.remainder(ntraces,2.0):
        print ('odd')
        npaires = int( (ntraces-1)/2)
    else:
        npaires = int(ntraces/2)
    
    # combine 10kV and 40kV data (assuming two adjacent columns)
    k=0
    merged_data = np.zeros((nsample,npaires))
    
    for i in np.arange(0,npaires*2,2):
        merged_data[:,k]=np.sum(traces[:,i:i+1],axis=1) 
        print(i,i+1)
        k = k+1
    
    
    
    for traceno in np.arange(0,npaires):
        # Prepare our data to be used in the filter #####################
        trace = merged_data[:,traceno] # get the proper trace in our file 
        ch = np.linspace(1,nsample,num=nsample) # Assign channel numbers 
        ch = ch[~np.isnan(trace)] # ignore the channels where there is no data
        trace = trace[~np.isnan(trace)] # ignore traces where there is no data
        np.savetxt(rdir+'CSV/merged/'+fname+'-paire-'+str(traceno)+'-merged-raw'+'.csv',np.transpose([ch,trace]),delimiter=',')
        
        ynoise, trace_clean = gigul.remove_background(ns,scale,trace,o,rdir + 'CSV/denoised/'+fname+'-paire-'+str(traceno)+'-denoised',ch)
        gigul.show_clean_trace (ch,trace,ynoise,trace_clean,rdir+'PNG/background/'+fname+'_paire'+str(traceno))
        plt.grid()
        plt.savefig(fname+'.png',format='png', transparent=True)     
        # # sommation des coups
        # T = 0  
        # trace_sum = np.zeros((1,npaires))
        
        # for T in merged_data[:,T]:
        #     T = 1
        #     trace_sum[:,T]= np.sum(merged_data[:,T],axis=-1)
        #     print (trace_sum)
        
        #     while trace_sum != 0:
        #        T = (T+1)
        #        trace_sum[:,T] = np.sum(merged_data[:,T:T+1],axis=-1)
        #        print (trace_sum)
        #     if trace_sum == 0:
        #        break
        # print ('stop')
    if a < b:
        continue
    if a == b:
        break
    
           
              
    
    #################################################################


    
plt.close('all')

