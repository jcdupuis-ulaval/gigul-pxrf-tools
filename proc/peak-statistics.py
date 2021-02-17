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
nbins = 128        # Number of bins that we want to establish  
imputation = -999  # If the peak is not found this is the value that gets added 
#########################################################################

data, FFID, n = gigul.merge_peak_data(pdir) # Merge the peak datasets into one and return the FFID and the number of datasets that were merged
trsh = int(n*0.6)  # Number of peaks that need to be detected in dataset to be relevant (peak should be present in 60% of the datasets)
gigul.generate_table_for_pca(data,fout,nbins,trsh,imputation,FFID) # Generate the PCA tables with peaks that were observed in more than 60% of the datasets