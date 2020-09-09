# -*- coding: utf-8 -*-
"""
Created on Wed Sep 09 2020
Routine to compute the statistics about the peaks that were found in the previous steps. This information will be used to 
bin the data and to see if some of the peaks are more prevalent than others 

@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# File setup for data ################################
pdir = '../results/peaks/'
# File list of all our directories ######################################
flist_peaks = os.listdir(path=pdir) # Peaks that were identified 
#########################################################################
i=0  # counter required to keep track of the matching files in our loop

for fname in flist_peaks:
    # We want to construct a database of peaks that were measured for the different measurements that were done 
    # Let's look at the files that are present 
    print ('Retrieving data from file :  %s' %fname)
    peaks = np.genfromtxt(pdir+flist_peaks[i],delimiter=',')
    i=i+1
