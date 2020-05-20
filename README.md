# Gigul_pxrf_tools
The tools in this module allows the user to process portable x-ray 
fluorescence data (pXRF) in order to identify elemental peaks. The data files
used as input are in CSV format. The first column contains the channel numbers
while the following columns contain the recorded spectra for each reading. It
is possible to process individual spectra or to merge adjacent spectra if they
were, for instance, acquired at different energy levels. 

These codes are still under early developpement and we make no representations
that they are fit for purpose. Use at your own risk.  
 
## Directory structure
The code, located in the *proc* folder, assumes a given folder structure. You 
may be required to set it up manually if your OS permissions do not allow 
the creation of new folders from your Python environment. The *data* folder
contains the raw datasets. The processing results will be stored in the *results*
folder. This folder contains the *CSV* folder, where the numerical processing
results are found and the *PNG* folder where the graphical representations are
found. The *peaks* folder contains the peaks that were estimated after processing. 

## Usual processing flow
Sample code helps the user get acquainted with the Gigul_pxrf_tools module. The first script to
run is the *remove-background-composite.py*. The name of the raw data file should be modified
by the user to match the name of the dataset to be processed. The *remove-background-composite.py*
merges adjacent columns for data that were acquired at two different energy levels. If your data was
acquired at only one energy, you can use *remove-background-individual-traces.py*. The results will be
stored in *../CSV/denoised/* and *../CSV/merged/*. The user can then run *pick-peaks.py* to identify 
the spectral peaks on the denoised data. The usermust select appropriate *amplitude and slope thersholds*
and the *peak width* of the peaks that are to be fitted. Figures will appear for each
of the files that are found in the *../CSV/denoised/* folder. The red lines that
appear on the raw spectra and the 1st derivative data allows to visualize the thresholds selected for the 
peak search algorithm that is used. The strategy implemented is adapted from "A Pragmatic Introduction to Signal Processing",
created and maintained by Prof. Tom O'Haver (https://terpconnect.umd.edu/~toh/spectrum/Differentiation.html)