## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved
## 24/11/2014

# IN : full tox dataset
# OUT : predictions for each compound



### Initialisations ##########

# Libraries
library(plyr)

#Path to the project directory
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_SCRIPTS_DIR='/home/az338/scripts/'

# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
TOX='BW_decrs'

# number of features
nF = 30

# number of down-sampling of the negative classification
# before classification
nDS = 100

# Load full dataset
load(file.path(CALC_DATA_DIR,"/datasets/fullDataSet_BW_decrs.Rdata"))

#  load leave_one_out script
source(file.path(CALC_SCRIPTS_DIR,"leave_one_out.R"))

##### ANALYSIS ############

## Random Forest
rf.predC = try(lv1out(fullDataset,nF,nDS,algo='rf',BioAc=F,GE=F)) # fingerprints only
rf.predBioA = try(lv1out(fullDataset,nF,nDS,algo='rf',FP=F,GE=F)) # bioactivity data only
rf.predG = try(lv1out(fullDataset,nF,nDS,algo='rf',FP=F,BioAc=F)) #gene expression only
rf.predCBioA = try(lv1out(fullDataset,nF,nDS,algo='rf',GE=F)) #fingerprints + bioactivity
rf.predGBioA = try(lv1out(fullDataset,nF,nDS,algo='rf',FP=F)) #gene expression + bioactivity
rf.predCG = try(lv1out(fullDataset,nF,nDS,algo='rf',BioAc=F)) #fingerprints + gene expression
rf.pred = lv1out(fullDataset,nF,nDS,algo='rf')  #all data sources

# save results
save(rf.pred,file="data/pred/rf.pred.Rdata")
save(rf.predC,file="data/pred/rf.predC.Rdata")
save(rf.predG,file="data/pred/rf.predG.Rdata")
save(rf.predBioA,file="data/pred/rf.predBioA.Rdata")
save(rf.predCBioA,file="data/pred/rf.predCBioA.Rdata")
save(rf.predGBioA,file="data/pred/rf.predGBioA.Rdata")
save(rf.predCG,file="data/pred/rf.predCG.Rdata")


## SVM
svm.predC = try(lv1out(fullDataset,nF,nDS,algo='svm',BioAc=F,GE=F)) # fingerprints only
svm.predBioA = try(lv1out(fullDataset,nF,nDS,algo='svm',FP=F,GE=F)) # bioactivity data only
svm.predG = try(lv1out(fullDataset,nF,nDS,algo='svm',FP=F,BioAc=F)) #gene expression only
svm.predCBioA = try(lv1out(fullDataset,nF,nDS,algo='svm',GE=F)) #fingerprints + bioactivity
svm.predGBioA = try(lv1out(fullDataset,nF,nDS,algo='svm',FP=F)) #gene expression + bioactivity
svm.predCG = try(lv1out(fullDataset,nF,nDS,algo='svm',BioAc=F)) #fingerprints + gene expression
svm.pred = lv1out(fullDataset,nF,nDS,algo='svm')  #all data sources

# save results
save(svm.pred,file="data/pred/svm.pred.Rdata")
save(svm.predC,file="data/pred/svm.predC.Rdata")
save(svm.predG,file="data/pred/svm.predG.Rdata")
save(svm.predBioA,file="data/pred/svm.predBioA.Rdata")
save(svm.predCBioA,file="data/pred/svm.predCBioA.Rdata")
save(svm.predGBioA,file="data/pred/svm.predGBioA.Rdata")
save(svm.predCG,file="data/pred/svm.predCG.Rdata")


