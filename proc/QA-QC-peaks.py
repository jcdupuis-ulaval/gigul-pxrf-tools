# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:49:33 2020
Routine to QA-QC the picks relative to the original dataset
@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# File setup for data ################################
rdir = '../results/CSV/merged/'
pdir = '../results/peaks/'
cdir = '../results/CSV/denoised/'
idir = '../results/PNG/'
# File list of all our directories ######################################
flist_raw = os.listdir(path=rdir) # Merged raw traces 
flist_clean = os.listdir(path=cdir) # Merged traces without background
flist_peaks = os.listdir(path=pdir) # Peaks that were identified 
#########################################################################
i=0  # counter required to keep track of the matching files in our loop

for fname in flist_raw:
    plt.figure()
    raw_data=np.genfromtxt(rdir+fname,delimiter=',')
    peaks = np.genfromtxt(pdir+flist_peaks[i],delimiter=',')
    clean_data = np.genfromtxt(cdir+flist_clean[i],delimiter=',')
    plt.plot(raw_data[:,0],raw_data[:,1],peaks[:,0],peaks[:,1],'+',clean_data[:,0],clean_data[:,1])
    plt.grid()
    plt.show()
    i=i+1

