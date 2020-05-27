import numpy as np
import matplotlib.pyplot as plt
import gigul_pxrf_tools as gigul
cdir = "../results/CSV/denoised/"
fname = "denoised-TDCAEXD332727z-paire-0.csv"
print (cdir+fname)
clean_data = np.genfromtxt(cdir+fname,delimiter=',')
ch = clean_data[:,0]
trace = clean_data[:,1]
n = 600
scale = 0.1
o = 2

scaled_trace = gigul.scale_trace(trace)
background = gigul.estimate_background(n,scale,scaled_trace,o,ch)
threshold = np.ones(len(ch))*scale*np.std(scaled_trace)

plt.plot(ch,threshold,ch,background,ch,gigul.smooth(gigul.smooth(scaled_trace)))

