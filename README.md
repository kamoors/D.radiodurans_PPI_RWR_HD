# D.radiodurans_PPI_RWR_HD
This repository contains the R code and Cytoscape files that were used for the project in the MSB course Network Biology.

Project code.R contains the R script that was used to extract the giant component from the filtered networks.

The following 3 files contain the full subnetworks for each of the proteins (integer in title correspons to confidence cut-off; e.g. 900 = 0.9 and above)

Analysis of RecA - 950_Deinococcus_radiodurans.sif
Analysis of UvrC - 900_Deinococcus_radiodurans.sif
Analysis of PprA - 700_Deinococcus_radiodurans.sif


The following 3 files contain the giant component for each of the subnetworks. The Random walk with restart and Heat diffusion algorithms were performed using these subnetworks:

Analysis of RecA - N950_giant.sif
Analysis of UvrC - N900_giant.sif
Analysis of PprA - N700_giant.sif


All plots can be generated by running the Project code.R file

All network visualizations are also available as files in the repository
