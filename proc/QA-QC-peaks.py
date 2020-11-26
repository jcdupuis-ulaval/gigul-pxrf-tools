# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:49:33 2020
Routine to QA-QC the picks relative to the original dataset
@author: chdup58
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# File setup for data and results################################
tdir = '../data/Terrain/'
rdir = '../results/CSV/merged/'
pdir = '../results/peaks/'
cdir = '../results/CSV/denoised/'
idir = '../results/PNG/QA_QC/'
fname_list = os.listdir(path=tdir)
b = len(fname_list)
a = 0
for value in fname_list:
    fname= fname_list[a]
    fname=fname[:-4]
    a = a+1
# Filter parameters #############################################
    flist_raw = os.listdir(path=rdir)
    flist_raw = [i for i in flist_raw if i.startswith(fname)]
    flist_clean = os.listdir(path=cdir)
    flist_clean = [i for i in flist_clean if i.startswith(fname)]
    flist_peaks = os.listdir(path=pdir)
    flist_peaks = [i for i in flist_peaks if i.startswith('picks_'+fname)]
    i=0
    for fname in flist_raw:
        plt.figure()
        raw_data=np.genfromtxt(rdir+fname,delimiter=',')
        peaks = np.genfromtxt(pdir+flist_peaks[i],delimiter=',')
        clean_data = np.genfromtxt(cdir+flist_clean[i],delimiter=',')
        plt.plot(raw_data[:,0],raw_data[:,1],peaks[:,0],peaks[:,1],'+',clean_data[:,0],clean_data[:,1])
        plt.grid()
        plt.title(fname)
        plt.legend(['original signal','peaks','normalised signal'])
        plt.xlabel('Channels')
        plt.ylabel('Total Count')
        plt.savefig(idir+fname[:-3]+'.png', format='png')
        plt.show()
        i=i+1
    if a < b:
        continue
    if a == b:
        break
