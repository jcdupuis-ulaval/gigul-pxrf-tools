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
ns=90                # Width of the window for noise estimate
scale = 0.05        # SNR Threshold
o = 1               # Order of the noise approximation 
#################################################################
# load the data file 
print ('Processing file : '+ddir+fname+'.csv')
data=np.genfromtxt(ddir+fname+'.csv',delimiter=',',skip_header=1)
m,n = data.shape

# Merge data
merged_data = np.sum(data[:,1:n],axis=1)


# first column is expected to be energy bins, start at second
'''
for traceno in np.arange(1,n):
    # Prepare our data to be used in the filter #####################
    trace = data[:,traceno] # get the proper trace in our file 
    ch = np.linspace(0,len(trace),num=len(trace)) # Assign channel numbers 
    ch = ch[~np.isnan(trace)] # ignore the channels where there is no data
    trace = trace[~np.isnan(trace)] # ignore traces where there is no data
#################################################################

    ynoise, trace_clean = gigul.remove_background(ns,scale,trace,o,rdir + 'CSV/'+fname+'_tr'+str(traceno),ch)
    gigul.show_clean_trace (ch,trace,ynoise,trace_clean,rdir+'PNG/'+fname+'_tr'+str(traceno))
    
plt.close('all')
'''