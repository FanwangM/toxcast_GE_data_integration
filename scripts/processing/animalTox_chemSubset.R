## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 6/11/14
## INPUT: toxcast animatl tox endpoint matrix
## OUTPUT: number of compounds from each dataset


CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'



library(plyr)
library(xtable)
library(gdata)


dat = read.csv(file.path(CALC_DATA_DIR,"toxcast_2014/animal_toxtoxrefdb_endpoint_matrix_AUG2014_FOR_PUBLIC_RELEASE.csv"),h=T)
chemInfo  = strsplit(as.character(dat[,1]),"\\|")
chemDtset = unlist(llply(chemInfo, function(c) return(unlist(c)[4])))
table(chemDtset)
