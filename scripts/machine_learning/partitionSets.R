
## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
#23/11/2014
# IN : full dataset (fingerprints + gene exp + bioact + tox outcome)
# OUT : training + test sets 

### Initialisations ##########

#Path to the project directory 
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


# Libraries
library(plyr)
library(gplots)
library(reshape2)

# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
TOX='BW_decrs'

# Full dataset
load(file.path(CALC_DATA_DIR,"datasets",paste("fullDataSet_",TOX,".Rdata",sep="")))

### Generate training set/ testing set ##########

# get positives (toxic)
tox.cp = which(fullDataset$liverCarcino == 1)
ntox.cp = which(fullDataset$liverCarcino == 0)

# Sample 1/4 of compounds with toxic outcoume
test.tox = sample(tox.cp,as.integer(length(tox.cp)/4))
# Sample 1/4 of compounds with non-toxic outcoume
test.ntox = sample(ntox.cp,as.integer(length(ntox.cp)/4))

# Test set = Compounds above
# sample is used here to shuffle the data
testSet =  fullDataset[sample(c(test.tox,test.ntox)),]

# Training set = All the compounds not selected above
trainingSet = fullDataset[-sample(c(test.tox,test.ntox)),]

# Save both training and test datasets
save(trainingSet,file=file.path(CALC_DATA_DIR,"datasets",paste("trainingSet_",TOX,".Rdata",sep="")))
save(testSet,file=file.path(CALC_DATA_DIR,"datasets",paste("testSet_",TOX,".Rdata",sep="")))

