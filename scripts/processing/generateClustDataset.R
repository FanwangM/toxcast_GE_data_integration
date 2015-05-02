
## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 28/11/14
## INPUT: activity matrix + in-vivo tox outcome + gene expression matrix 
## OUTPUT: training set and testing set 


### Initialisations ##########

#Path to the project directory 
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


# Libraries
library(plyr)
library(reshape2)
library(impute)

#Cell line for gene expression dataset
CLINE = "MCF7"

# load boolean activity matrix (bioactivity target)
load(file.path(CALC_DATA_DIR,"datasets/full_boolAct_mat.Rdata"))
load(file.path(CALC_DATA_DIR,"datasets/full_boolAct_anno.Rdata"))

## load fingerprints
fp_df = read.table(file.path(CALC_DATA_DIR,"/datasets/fingerprints.tab"),h=F,sep="\t",quote="",comment.char="",row.names=1)
colnames(fp_df) = unlist(llply(1:ncol(fp_df), function(i) paste("FP_",i,sep="")))

## load gene expression dataset (name:[CLINE]_df,type:df)
load(file.path(CALC_DATA_DIR,"boolGE",paste("boolean_geneExp_",CLINE,".Rdata",sep="")))
GE_df = get(ls()[grep(paste(CLINE,"_df",sep=""),ls())])

## load animalTox matrix (boolean activity)
load(file.path(CALC_DATA_DIR,"datasets/bool_animalTox.Rdata"))
# contingency table of chemicals and effects
cttb = ddply(animaltox,.(chid,effect),summarize,n=length(effect))
cttb = acast(cttb,chid~effect,value.var="n")
cttb[which(is.na(cttb))] = 0
cttb[which(cttb > 0)] = 1

# remove column corresponding to no effect(?)
animaltox_mat= cttb[,-1]

# load chemicals dataset to map the chemicals
# in the above two datasets ('linker' dataset)
load(file.path(CALC_DATA_DIR,"chemicals.Rdata"))

### Dataset Fusion ##########

# match animal tox endpoints with bioactivity
tox.ind =  match(rownames(animaltox_mat),rownames(boolAct_mat))
boolAct_mat = boolAct_mat[tox.ind,]
boolAct_mat = data.frame(boolAct_mat,animaltox_mat)

# match toxcast ids in the activity matrix with those of the chemical dataset
cmn_chems = match(rownames(boolAct_mat),chemicals$TC_ID)

# reduce gene expression dataset to chemicals
# represented in the tox dataset
GE_df = subset(GE_df,rownames(GE_df) %in% chemicals$LINCS_ID[cmn_chems])

# reorder gene expression dataset to match order of tox dataset
GE_df = GE_df[match(chemicals$LINCS_ID[cmn_chems],rownames(GE_df)),]

# reorder fingerprint dataset to match order of gene expression dataset
fp_df = fp_df[match(rownames(GE_df),rownames(fp_df)),]

# clustering dataset
clusDataset = data.frame(boolAct_mat,GE_df,fp_df)

# impute missing values
# Store position of missing values 
nas = which(is.na(clusDataset),arr.ind=T)

# Impute missing values using KNN
clusDataset = impute.knn(as.matrix(clusDataset))[[1]]

# Because KNN impute via average, results are transformed
# back into binairy variables by rounding the values
clusDataset[nas] = round(clusDataset[nas])


save(clusDataset,file=file.path(CALC_DATA_DIR,"datasets/clustDataset.Rdata"))
