## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved
## 8/12/2014

### Initialisations ##########

#Path to the project directory
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_SCRIPTS_DIR='/home/az338/scripts/'


# Libraries
library(plyr)


# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
#TOX = "livrCarcino"
TOX='BW_decrs'


# number of down-sampling of the negative classification
# before classification
nDS = 100

# Load full dataset (name: fullDataset, type:df)
load(file.path(CALC_DATA_DIR,"/datasets/fullDataSet_BW_decrs.Rdata"))

#  load leave_one_out script (function: lv1out)
source(file.path(CALC_SCRIPTS_DIR,"leave_one_out.R"))

# Number of features tested
nF = c(5,10,20,50)

# Combinations of data tested
DataCombi = list(
    FP=c(T,F,F),
    BioAc=c(F,T,F),
    GE=c(F,F,T),
    FP.BioAc=c(T,T,F),
    FP.GE=c(T,F,T),
    BioAc.GE=c(F,T,T),
    ALL=c(T,T,T))

# Loop over combinations of data
rfRes=llply(1:length(DataCombi), function(i) {
    cat("Combination:", names(DataCombi[i]),'\n')
    # loop over numer of features
    return(llply(nF, function(nF) {
        cat("Number of features tested:",nF,'\n')
        # and get random forest predictions
        return(lv1out(fullDataset,nF,nDS,'rf',FP=DataCombi[[i]][1],BioAc=DataCombi[[i]][2],GE=DataCombi[[i]][3]))
    }))
})
save(rfRes,file=file.path(CALC_DATA_DIR,'pred/RF_fs_DCombi_preds.Rdata'))

# Loop over combinations of data
svmRes=llply(1:length(DataCombi), function(i) {
    cat("Combination:", names(DataCombi[i]),'\n')
    # loop over numer of features
    return(llply(nF, function(nF) {
        cat("Number of features tested:",nF,'\n')
        # and get svm predictions
        return(lv1out(fullDataset,nF,nDS,'svm',FP=DataCombi[[i]][1],BioAc=DataCombi[[i]][2],GE=DataCombi[[i]][3]))
    }))
})
save(svmRes,file=file.path(CALC_DATA_DIR,'pred/SVM_fs_DCombi_preds.Rdata'))



