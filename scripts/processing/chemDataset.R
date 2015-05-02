## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 17/10/14
## INPUT: LINCS+TOXCAST InchiKeys
## OUTPUT: Chemical dataset with compounds shared by LINCS/ToxCast


PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the data for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


library(plyr)
library(xtable)
library(gdata)



# DSSTOX imports
# Load standardized compounds ids for dsstox
# and corresponding smiles
tb=read.table(
        file.path(DATA_DIR,"toxcast_2013/chemicals/DSSTox_All_Smiles_Standardized.tab"),h=F,stringsAsFactors=F,row.names=1)
dsstox= tb[,c(4,1)]
rm(tb)
# Load correspondong InchiKeys
tb=read.table(file.path(DATA_DIR,"toxcast_2013/chemicals/DSSTox_InchiKeys.tab"),
                    sep='\t',h=F,stringsAsFactors = F,row.names=1)
tb$Inchikey=unlist(llply(tb[,1],
                            function(x) unlist(strsplit(unlist(strsplit(x,"\n"))[2],
                                                        "InChIKey="))[2]))
dsstox$inchikey = tb$Inchikey
colnames(dsstox) = c("GSID","SMILES","INCHIKEY")
rm(tb) 

# TOXCAST imports
# Load standardized compounds ids for toxcast
# and corresponding SMILES
toxcast=read.table(
        file.path(DATA_DIR,"toxcast_2013/chemicals/ToxCast_All_Smiles_Standardized.tab"),h=F,sep="\t",stringsAsFactors=F,row.names=1)
toxcast = toxcast[,c(4,1)]
# Load corresponding Inchikeys
tb = read.table(file.path(DATA_DIR,"toxcast_2013/chemicals/ToxCast_InchiKeys.tab"),
                     sep='\t',h=F,stringsAsFactors = F,row.names=1)
tb$Inchikey=unlist(llply(tb[,1],
                            function(x) unlist(strsplit(unlist(strsplit(x,"\n"))[2],
                                                        "InChIKey="))[2]))
toxcast$inchikey = tb$Inchikey
colnames(toxcast) = c("GSID","SMILES","INCHIKEY")
rm(tb)

#LINCS imports
# Load standardized compounds ids for LINCS
# and corresponding SMILES
lincs = read.table(
        file.path(DATA_DIR,"lincs/All_Compounds_Smiles_Standardized.tab"),h=F,sep="\t",stringsAsFactors=F,row.names=1)
lincs = lincs[,c(4,1)]
tb =  read.table(file.path(DATA_DIR,"lincs/All_Compounds_InchiKeys.tab"),
                     sep='\t',h=F,stringsAsFactors = F,row.names=1)
tb$Inchikey=unlist(llply(tb[,1],
                            function(x) unlist(strsplit(unlist(strsplit(x,"\n"))[2],
                                                        "InChIKey="))[2]))
lincs$INCHIKEY = tb$Inchikey
colnames(lincs) = c("LINCS_ID","SMILES","INCHIKEY")


# Remove duplicated compounds
lincs = unique(lincs)
toxcast = unique(toxcast)
dsstox = unique(dsstox)
dsstox=dsstox[-which(dsstox$INCHIKEY %in% toxcast$INCHIKEY),]

#The list of shared Inchikeys (obtained in Overlap\_Analysis.R) was used
#to link the LINCS compounds to the ToxCast/DSSTox datasets.


# Load shared InchiKeys 
cmn_Inchikeys = readLines(file.path(DATA_DIR,"intermediate_files/cmn_Inchikeys"))

# Eliminate compounds which 
#do not have an inchikey in that list
dsstox=subset(dsstox,INCHIKEY %in% cmn_Inchikeys)
toxcast=subset(toxcast,INCHIKEY %in% cmn_Inchikeys)
lincs=subset(lincs,INCHIKEY%in% cmn_Inchikeys)

# Match LINCS ids to ToxCast ids
lincs$TC_ID=toxcast[match(lincs$INCHIKEY,toxcast$INCHIKEY),"GSID"]

# Match LINCS ids to DSSTox ids
ids = dsstox[match(lincs$INCHIKEY,dsstox$INCHIKEY),"GSID"]
# Reformat the gsid to match the format in the primary toxcast dataset
lincs$TC_ID[which(!is.na(ids))] = paste("DSSTox_GSID_",ids[which(!is.na(ids))],sep="")


#Export to Rdata object
chemicals = lincs[,]
save(chemicals,file=file.path(DATA_DIR,"/intermediate_files/chemicals.Rdata"))

# save a copy to calculon as well
system(paste("scp ",file.path(DATA_DIR,"/intermediate_files/chemicals.Rdata")," calculon:",file.path(CALC_DATA_DIR,"boolGE"),sep=""))
