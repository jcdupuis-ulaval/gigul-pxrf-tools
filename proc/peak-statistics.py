# -*- coding: utf-8 -*-
"""
Created on Wed Sep 09 2020
Routine to collate the peaks that were found in the previous steps and compute basic statistics on the channel number and 
the mean amplitude that were found previously. When peaks are repeatable (e.g. found in multiple measurements) they are saved
to a csv file  (fout) for later PCA analysis  

@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
import gigul_pxrf_tools as gigul
import csv
import os

# File setup for data ################################
pdir = '../results/peaks/'
# File list of all our directories ######################################
flist_peaks = os.listdir(path=pdir) # Peaks that were identified 
fout = pdir+'peaks-for-PCA'
#########################################################################
i=0  # counter required to keep track of the matching files in our loop
ch = []
amp = []
fid = []
trsh = 4

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

# Compute the stastics on the peaks that we have identified
mu_ch,std_ch,mu_amp,std_amp,n = gigul.calc_stats(data)
# Compare the average peaks to the peaks that were identified 
plt.subplot(2,1,1)
plt.plot(mu_ch,mu_amp,'+')
plt.plot(data[:,1],data[:,2],'.')
plt.legend('mean','originals')
plt.subplot(2,1,2)
plt.plot(mu_ch,n,'.')
plt.ylabel('Number of points in mean')
plt.show()
print (n)

# For our PCA analysis we will want only the peaks that can be seen in all our datasets (e.g. are repeatable)

data_clean = np.hstack((mu_ch[np.where(n>trsh)[0]],
std_ch[np.where(n>trsh)[0]],
mu_amp[np.where(n>trsh)[0]],
std_amp[np.where(n>trsh)[0]]))

# Saves the data that is ready for PCA 
# C1 : Mean channel 
# C2 : STD channel
# C3 : Mean Amplitude
# C4 : STD Amplitude
np.savetxt(fout+'.csv',data_clean,delimiter=',')

