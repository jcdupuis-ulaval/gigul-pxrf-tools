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

# The channel numbers may not be integers because of the picking process, we reduce the data to integer channels only 
ch = np.round(np.array(ch).astype(float),decimals=0)
amp = np.array(amp).astype(float)
fid = np.array(fid).astype(int)
values = np.hstack((np.vstack(fid),np.vstack(ch),np.vstack(amp)))
dtype = [('fid',int),('ch',float),('amp',float)]
data = np.array(values,dtype=dtype)

data_sorted = np.sort(data,order='ch')
m,n = data_sorted.shape
for i in np.arange(m):
    print (i)
#df = pd.DataFrame(data,columns = ['fid','ch','amp'])

print (data)

# At this stage we have all of the peaks that were found 
plt.scatter(ch,amp,c=fid)
plt.show()

