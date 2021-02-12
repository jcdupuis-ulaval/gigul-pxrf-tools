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
odir = '../results/PCA-tables/'
# File list of all our directories ######################################
flist_peaks = os.listdir(path=pdir) # Peaks that were identified 
fout = odir+'PCA-mean.csv'
#########################################################################
trsh = 4
nbins = 128
imputation = -999
#########################################################################
i=0  # counter required to keep track of the matching files in our loop
ch = []
amp = []
fid = []
FFID = []
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
        # We want to save the name of the original file and the pair of readings 
        ffid = fname.split("-")
        FFID.append(ffid[1] +'-' + str(i))
        i=i+1

ch = np.round(np.array(ch).astype(float),decimals=2)
amp = np.array(amp).astype(float)
fid = np.array(fid).astype(int)

# At this stage we have all of the peaks that were found and stored in files 
#fig, ax = plt.subplots()
#scatter = ax.scatter(ch,amp,c=fid)
#legend1 = ax.legend(*scatter.legend_elements(),title='Paires')
#ax.add_artist(legend1)

#plt.show()

# Combine the data in one array 
data = np.hstack((np.vstack(fid),np.vstack(ch),np.vstack(amp)))
#dtype = [('fid',int),('ch',float),('amp',float)]
#data = np.array(values,dtype=dtype)

# Now that we have our data we can count the number of occurences in each bin 
# this can simply be done using plt.hist

x = plt.hist(data[:,1],nbins)

# We want to find the bins where the number of observations exceeds our threshold value 
idx = np.where(x[0]>trsh)
m,n = np.shape(idx)


# Now that we have the indices of the bins that exceed the thershold we want to split our data back into their individual
# pairs of measurements and determine the average amplitudes for the bins where we expect to have data according to the bin edges

#pca_file = open(fout,"w")
pca_matrix_mean = np.zeros(((len(FFID),n-1)))
pca_matrix_std = np.zeros(((len(FFID),n-1)))
pca_matrix_bins = np.zeros(((len(FFID),n-1)))
for i in np.arange(len(FFID)):
    meas = data[np.where(data[:,0]==i)]
    print ('Looking into dataset %s' %FFID[i])
    for j in np.arange(n-1):
        lower_edge = x[1][idx[0][j]]
        upper_edge = x[1][idx[0][j+1]]
        pca_matrix_bins [i,j]= lower_edge + (upper_edge-lower_edge)/2.0
        # print our bin edges so that we can see QA where they fall
        print('Selecting from %2.4f to %2.4f' %(lower_edge,upper_edge))
        amp = meas[np.where((meas[:,1]>=lower_edge) & (meas[:,1]<=upper_edge )),1]
        if amp.size>1: # we have many points that fall in this bin
            pca_matrix_mean[i,j] = np.mean(meas[np.where((meas[:,1]>=lower_edge) & (meas[:,1]<=upper_edge )),2]) # compute the average amplitude of the points that fall in this bin
            pca_matrix_std[i,j] = np.std(meas[np.where((meas[:,1]>=lower_edge) & (meas[:,1]<=upper_edge )),2])  # compute the std of the points that fall in this bin 
        elif amp.size == 0.0: # We did not find a peak that has this value in this dataset - enter imputation value 
            pca_matrix_mean[i,j] = imputation
        else:
            pca_matrix_mean[i,j] = meas[np.where((meas[:,1]>=lower_edge) & (meas[:,1]<=upper_edge )),2]
            pca_matrix_std[i,j] = np.nan
            
        print('found this peak %2.4f' %pca_matrix_bins[i,j])
        print('amplitude is  %2.4f' %pca_matrix_mean[i,j])
        print ('Saving in bin %2.4f' %pca_matrix_bins[i,j])
        #file.write (%s,#2.4f,)


# Now that we have all of the data in our matrices we want to write these to files that can be used for PCA analysis
# Since PCA could be done with other software, it's useful to put a header at the top of the file and identify each measurement
# with the analysis number. 

f = open(fout,'w')

hdr = np.array2string(pca_matrix_bins[0,:],separator=',')
hdr = hdr.replace('[','Reading-ID,')
hdr = hdr.replace(']','')
hdr = hdr.replace('\n','')
f.write(hdr+'\n')

for i in np.arange(len(FFID)):

    reading = np.array2string(pca_matrix_mean[i,:],separator=',')
    reading = reading.replace('[',FFID[i]+',')
    reading = reading.replace('\n','')
    reading = reading.replace(']','')
    f.write(reading+'\n')

f.close()

print ()
'''
pca_file = open(fout,'w')
pca_file.write(pca_matrix_bins[0,:].tostring())
pca_file.write(pca_matrix_mean[0,:].tostring())
pca_file.close()
print (pca_matrix_mean)
'''