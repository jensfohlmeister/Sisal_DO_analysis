# Sisal_DO_analysis
This code belongs to the following paper: Fohlmeister et al., 2023, PNAS

################################################################################
# this program identifies D/O events of some predefined records from SISAL2
# The identified D/O events are then stacked and the mean offset of stadial and
# interstadial conditions is provided as output.
#
# This work is intended to be published in PNAS
#
# To run this model, please download and unzip the SISAL2 data base from here:
# https://researchdata.reading.ac.uk/256/
# and save the SISAL2-files in a './SISALv2_csv' file folder 
# the best DO bearing speleothem data covering the last glacial were prescreened and 
# are provided in DO-direction.txt, which also includes information about the 
# direction of stable oxygen isotope shifts based on information of the individual
# original data papers.
#
# Furthermore, the program will need NGRIP data (as published in Rasmussen et 
# al., 2014), which are also provided in two files in this file folder for 
# your convenience (the files should be located in the same directory as the 
# julia file)
################################################################################
