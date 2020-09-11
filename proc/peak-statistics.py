# -*- coding: utf-8 -*-
"""
Created on Wed Sep 09 2020
Routine to compute the statistics about the peaks that were found in the previous steps. This information will be used to 
bin the data and to see if some of the peaks are more prevalent than others 

@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os

# File setup for data ################################
pdir = '../results/peaks/'
# File list of all our directories ######################################
flist_peaks = os.listdir(path=pdir) # Peaks that were identified 
#########################################################################
i=0  # counter required to keep track of the matching files in our loop
ch = []
amp = []
fid = []

for fname in flist_peaks:
    # We want to construct a database of peaks that were measured for the different measurements that were done 
    # Let's look at the files that are present 
    print ('Retrieving data from file :  %s' %fname)
    with open(pdir+flist_peaks[i]) as data_in:
        csv_reader = csv.reader(data_in,delimiter=',')
        line_count = 0
        for row in csv_reader:
            fid.append(i)
            ch.append(row[0])
            amp.append(row[1])
        
        i=i+1

ch = np.round(np.array(ch).astype(float),decimals=2)
amp = np.array(amp).astype(float)
fid = np.array(fid).astype(int)

# At this stage we have all of the peaks that were found and stored in files 
fig, ax = plt.subplots()
scatter = ax.scatter(ch,amp,c=fid)
legend1 = ax.legend(*scatter.legend_elements(),title='Paires')
ax.add_artist(legend1)

plt.show()

# Combine the data in one array 
data = np.hstack((np.vstack(fid),np.vstack(ch),np.vstack(amp)))
#dtype = [('fid',int),('ch',float),('amp',float)]
#data = np.array(values,dtype=dtype)

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
        ch = data[ch_edges[i]+1:ch_edges[i+1],1]
        amp = data[ch_edges[i]+1:ch_edges[i+1],2]
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

mu_ch,std_ch,mu_amp,std_amp,n = calc_stats(data)
plt.subplot(2,1,1)
plt.plot(mu_ch,mu_amp,'+')
plt.plot(data[:,1],data[:,2],'.')
plt.subplot(2,1,2)
plt.plot(mu_ch,n,'.')
plt.show()
print (n)
'''
data_sorted = np.sort(data,order='ch')
m,n = data_sorted.shape
for i in np.arange(m):
    print (i)
#df = pd.DataFrame(data,columns = ['fid','ch','amp'])

print (data)
'''
# At this stage we have all of the peaks that were found 
#plt.scatter(ch,amp,c=fid)
#plt.show()

