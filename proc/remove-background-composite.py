# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:44:00 2020

@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
import gigul_pxrf_tools as gigul


# File setup for data and results################################
fname ='TDCAEXD332727z'
ddir = '../data/'
rdir = '../results/'

# Filter parameters #############################################
ns=50               # Width of the window for noise estimate
scale = 0.05        # SNR Threshold
o = 1               # Order of the noise approximation 
#################################################################
# load the data file 
print ('Processing file : '+ddir+fname+'.csv')
data=np.genfromtxt(ddir+fname+'.csv',delimiter=',',skip_header=1)
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
    
    np.savetxt(rdir+'CSV/merged/'+'merged-raw-'+fname+'-paire-'+str(traceno)+'.csv',np.transpose([ch,trace]),delimiter=',')
#################################################################

    ynoise, trace_clean = gigul.remove_background(ns,scale,trace,o,rdir + 'CSV/denoised/'+'denoised-'+fname+'-paire-'+str(traceno),ch)
    gigul.show_clean_trace (ch,trace,ynoise,trace_clean,rdir+'PNG/'+fname+'_paire'+str(traceno))
    
plt.close('all')
