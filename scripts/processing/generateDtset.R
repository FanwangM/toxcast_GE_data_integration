
## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 3/11/14
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
library(gplots)
library(reshape2)

# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
TOX='BW_decrs'

## load fingerprints
fp_df = read.table(file.path(CALC_DATA_DIR,"/datasets/fingerprints.tab"),h=F,sep="\t",quote="",comment.char="",row.names=1)
colnames(fp_df) = unlist(llply(1:ncol(fp_df), function(i) paste("FP_",i,sep="")))

## load gene expression dataset (name:[CLINE]_df,type:df)
load(file.path(CALC_DATA_DIR,"boolGE",paste("boolean_geneExp_",CLINE,".Rdata",sep="")))
GE_df = get(ls()[grep(paste(CLINE,"_df",sep=""),ls())])

# load bioactivity dataset with in-vivo tox outcome (name:boolAct_mat,type:df)
load(file.path(CALC_DATA_DIR,"datasets",paste("boolAct_mat_",TOX,".Rdata",sep="")))

# load chemicals dataset to map the chemicals
# in the above two datasets (linker dataset)
load(file.path(CALC_DATA_DIR,"chemicals.Rdata"))

### Dataset Fusion ##########

#Match toxcast ids in the activity matrix with those of the chemical dataset
cmn_chems = match(rownames(boolAct_mat),chemicals$TC_ID)

# reduce gene expression dataset to chemicals
# represented in the tox dataset
GE_df = subset(GE_df,rownames(GE_df) %in% chemicals$LINCS_ID[cmn_chems])

# reorder gene expression dataset to match order of tox dataset
GE_df = GE_df[match(chemicals$LINCS_ID[cmn_chems],rownames(GE_df)),]

# reorder fingerprint dataset to match order of gene expression dataset
fp_df = fp_df[match(rownames(GE_df),rownames(fp_df)),]

# add gene expression dataset to tox dataset
fullDataset = data.frame(boolAct_mat,GE_df,fp_df)

#heatmap.2(as.matrix(fullDataset),dendrogram="none",trace="none",col=c("blue","orange","red"))

# save dataset to R object
save(fullDataset,file=file.path(CALC_DATA_DIR,"datasets",paste("fullDataSet_",TOX,".Rdata",sep="")))
