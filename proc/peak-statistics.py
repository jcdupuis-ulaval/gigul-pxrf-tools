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

# File setup for data ##################################################
pdir = '../results/peaks/'
odir = '../results/PCA-tables/'
fout = odir+'PCA-mean.csv'
#########################################################################
trsh = 1           # Number of peaks that need to be detected in dataset to be relevant
nbins = 128        # Number of bins that we want to establish  
imputation = -999  # If the peak is not found this is the value that gets added 
#########################################################################

data, FFID = gigul.merge_peak_data(pdir)
gigul.generate_table_for_pca(data,fout,nbins,trsh,imputation,FFID)