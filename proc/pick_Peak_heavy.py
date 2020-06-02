import numpy as np
import gigul_pxrf_tools as gigul
import os

# File setup for data and results################################
#fname ='TDCAEXD332727z_paire5'
ddir = '../results/CSV/denoised/'
rdir = '../results/peaks/heavy/'
idir = '../results/PNG/heavy/'
# Filter parameters #############################################
amp_threshold=10.0
slope_threshold =-10.0
peak_half_width = 5
#################################################################

flist = os.listdir(path=ddir)
    
for fname in flist:
    data=np.genfromtxt(ddir+fname,delimiter=',')
    trace = data[1100:2048,1]
    ch = data[1100:2048:,0] 
    peak_est = gigul.estimate_peaks(trace,ch,amp_threshold,slope_threshold)
    print (peak_est)
    peaks = gigul.refine_peaks(trace,peak_est,peak_half_width,ch,rdir+'picks_'+fname)
    gigul.show_peaks(trace,ch,peaks,peak_est,idir+'picks_'+fname)
