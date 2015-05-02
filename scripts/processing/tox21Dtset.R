
## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 3/11/14
## INPUT: Tox21 dataset
## OUTPUT: gain-log(ac50) matrix + corresponding annotation of the points



### Initialisations ##########

#Path to the project directory 
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

TOX='BW_decrs'
TOX_GREP = "BodyWeight.Body Weight.Decrease" #toxicity query to match in tox matrix

# Libraries
library(plyr)
library(gplots)
library(reshape2)

# Load chemical dataset (name:chemicals,type:df)
load(file.path(CALC_DATA_DIR,"/chemicals.Rdata"))

# Read in Tox21 summary csv file
dat = read.csv(file.path(CALC_DATA_DIR,"toxcast_2014/ToxCast_Tox21_Level5&6_20141022.csv"),h=T)

# Reduce dataset to chemicals in the analysis (LINCS/ToxCast)
dat = subset(dat,chid %in% chemicals$TC_ID)

### Boolean activity matrix ##########

# Compute matrix of activity(1)/inactivity(0)
boolAct_mat = acast(dat,chid~aenm,value.var="hitc",fun.aggregate=max)
boolAct_mat[which(boolAct_mat == -Inf)] = NA
save(boolAct_mat,file=file.path(CALC_DATA_DIR,"datasets/full_boolAct_mat.Rdata"))

# Corresponding annotation(names, flags etc...)
boolAct_anno = dat[,c("chid","chnm","aeid","aenm","modl","modl_rmse","nconc","npts","nrep","actp","resp_unit","flag")]
save(boolAct_anno,file=file.path(CALC_DATA_DIR,"datasets/full_boolAct_anno.Rdata"))

### Matrix of gain-log(AC50) values ##########

# Reduce dataset to active data only (hitcall == 1)
dat = subset(dat,hitc == 1)

# Compute gain-log(AC50) matrix (rows: chemicals, cols:assays)
# 594 chemicals for 61 assays
# miss.values = 87%
glac50_mat = acast(dat,chid~aenm,value.var="modl_ga",fun.aggregate=median)

# save gain-log(ac50) matrix
save(glac50_mat,file = file.path(CALC_DATA_DIR,"datasets/full_glac50_mat.Rdata"))

# save corresponding annotation of the chemical/assay points (names, flags etc...)
glac50_anno = dat[,c("chid","chnm","aeid","aenm","modl","modl_rmse","nconc","npts","nrep","actp","resp_unit","flag")]
save(glac50_anno,file = file.path(CALC_DATA_DIR,"datasets/full_glac50_anno.Rdata"))

### Animal toxicity with boolean bioactivity  ##########

# Read in animal tox csv file
animaltox = read.csv(file.path(CALC_DATA_DIR,"toxcast_2014/animal_tox/toxrefdb_study_tg_effect_endpoint_AUG2014_FOR_PUBLIC_RELEASE.csv"),h=T)

# Reformat the chemical id to match those in tox21 dataset
animaltox$chid = gsub("^DSSTox_GSID_","",animaltox$chemical_id)

# Get the full list of compounds ids represented in this matrix
anTox_ids = unique(animaltox$chid)

# Match animal toxicity with tox21 compounds on the boolean matrix
mtchg = match(rownames(boolAct_mat),anTox_ids)
length(which(!is.na(mtchg)))
#96 compounds


# Reduce animal tox dataset to these compounds
animaltox = subset(animaltox,chid %in% rownames(boolAct_mat) & species %in% c("mouse","rat"))

# Concatenate effect type, description and direction
animaltox= transform(animaltox,effect=interaction(endpoint_target,effect_desc,direction,sep="."))
save(animaltox,file=file.path(CALC_DATA_DIR,"datasets/bool_animalTox.Rdata"))

# contingency table of chemicals and effects
cttb = ddply(animaltox,.(chid,effect),summarize,n=length(effect))
cttb = acast(cttb,chid~effect,value.var="n")
cttb[which(is.na(cttb))] = 0
cttb[which(cttb > 0)] = 1

# remove column corresponding to no effect(?)
cttb= cttb[,-1]

# get top 20 more represented in vivo tox 
tail(sort(colSums(cttb)),20)

# heatmap
#heatmap.2(cttb,trace="none",col=c("orange","red"))

# select one type of toxicity
# here : bodyweight decrease (any)
col_idx= grep(TOX_GREP,colnames(cttb))


# Select compound matching this tox endpoint
# BE AWARE OF DUPLICATED COMPOUNDS (and eliminate duplicates)
if(length(col_idx) > 1)
    tox.cp = rownames(cttb)[which(cttb[,col_idx] == 1,arr.ind=T)[,1]]
else
    tox.cp = rownames(cttb)[which(cttb[,col_idx] == 1)]

# reduce dataset to chemicals that appear in the animal tox dataset only
boolAct_mat = subset(boolAct_mat, rownames(boolAct_mat) %in% anTox_ids)

# add tox endpoint to the boolean activity matrix
boolAct_mat = as.data.frame(boolAct_mat)
boolAct_mat$toxClass = 0 
boolAct_mat[tox.cp,"toxClass"] = 1

# replace column name named "toxClass" by real tox name
# specified by the constant at the top
indx = grep("toxClass",colnames(boolAct_mat))
colnames(boolAct_mat)[indx] = TOX

save(boolAct_mat,file=file.path(CALC_DATA_DIR,"datasets",paste("boolAct_mat_",TOX,".Rdata",sep="")))

# heatmap
#heatmap.2(as.matrix(boolAct_mat[order(boolAct_mat$liverCarcino),]),dendrogram="none",Rowv=NULL,Colv=NULL,trace="none",col=c("grey","red"))



### Animal toxicity with ac50 matrix ##########

# Read in animal tox csv file
animaltox = read.csv(file.path(CALC_DATA_DIR,"toxcast_2014/animal_tox/toxrefdb_study_tg_effect_endpoint_AUG2014_FOR_PUBLIC_RELEASE.csv"),h=T)

# Reformat the chemical id to match those in tox21 dataset
animaltox$chid = gsub("^DSSTox_GSID_","",animaltox$chemical_id)


# Match animal toxicity with tox21 compounds on the activity matrix
mtchg = match(rownames(glac50_mat),unique(animaltox$chid))
length(which(!is.na(mtchg)))
# 83 compounds


# Reduce animal tox dataset to these compounds
animaltox = subset(animaltox,chid %in% rownames(glac50_mat) & species %in% c("mouse","rat"))

# Concatenate effect type, description and direction
animaltox= transform(animaltox,effect=interaction(endpoint_target,effect_desc,direction,sep="."))
save(animaltox,file=file.path(CALC_DATA_DIR,"datasets/glac50_animaltox.Rdata"))

# contingency table of chemicals and effects
cttb = ddply(animaltox,.(chid,effect),summarize,n=length(effect))
cttb = acast(cttb,chid~effect,value.var="n")
cttb[which(is.na(cttb))] = 0
cttb[which(cttb > 0)] = 1

# remove column corresponding to no effect(?)
cttb= cttb[,-1]


# get top 20 more represented in vivo tox 
tail(sort(colSums(cttb)),20)


# select one type of toxicity
# here : bodyweight decrease (any)
col_idx= grep(TOX_GREP,colnames(cttb))


# Select compound matching this tox endpoint
# BE AWARE OF DUPLICATED COMPOUNDS (and eliminate duplicates)
if(length(col_idx) > 1)
    tox.cp = rownames(cttb)[which(cttb[,col_idx] == 1,arr.ind=T)[,1]]
else
    tox.cp = rownames(cttb)[which(cttb[,col_idx] == 1)]



# reduce dataset to chemicals that appear in the animal tox dataset only
glac50_mat = subset(glac50_mat, rownames(glac50_mat) %in% anTox_ids)


glac50_mat = as.data.frame(glac50_mat)
glac50_mat$toxClass = 0 
glac50_mat[tox.cp,"toxClass"] = 1

# replace column name named "toxClass" by real tox name
# specified by the constant at the top
indx = grep("toxClass",colnames(glac50_mat))
colnames(glac50_mat)[indx] = TOX

save(glac50_mat,file=file.path(CALC_DATA_DIR,"datasets",paste("glac50_mat_",TOX,".Rdata",sep="")))


#heatmap
#heatmap.2(as.matrix(glac50_mat[order(glac50_mat$liverCarcino),]),dendrogram="none",Rowv=NULL,Colv=NULL,trace="none")

