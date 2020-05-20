# Gigul_pxrf_tools
The tools in this package allows the user to process portable x-ray 
fluorescence data (pXRF) in order to identify elemental peaks. The data files
used as input are in CSV format. The first column contains the channel numbers
while the following columns contain the recorded spectra for each reading. It
is possible to process individual spectra or to merge adjacent spectra if they
were, for instance, acquired at different energy levels. 

These codes are still under early developpement and we make no representations
that they are fit for purpose. Use at your own risk.  
 


Les données en format CSV sont dans le dossier : data
Les codes pour l'analyse sont dans le dossier : proc
Le script remove-background permet de réduire les effets des éléments légers
Le script pick-peaks permet d'identifier les événements
Les résultas de l'analyse se trouve dans : results
