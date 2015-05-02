## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 2/11/14
## INPUT: LINCS+TOXCAST(2014) InchiKeys
## OUTPUT: Chemical dataset with compounds shared by LINCS/ToxCast for the 2014 realease


PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the data for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

## It also executes libraries that are usually used by the
## scripts in the project

library(plyr)
library(xtable)
library(gdata)



# TOXCAST imports
# Load standardized compounds ids for toxcast
# and corresponding SMILES
toxcast=read.table(
        file.path(DATA_DIR,"toxcast_2014/chemicals/ToxCast_All_Smiles_Standardized/chemicals/ToxCast_SMILES2014_Standardized.tab"),h=F,sep="\t",stringsAsFactors=F,row.names=1)
toxcast = toxcast[,c(4,1)]
# Load corresponding Inchikeys
tb = read.table(file.path(DATA_DIR,"toxcast_2014/chemicals/ToxCast_InchiKeys.tab"),
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

# Match LINCS ids to ToxCast ids
lincs$TC_ID=toxcast[match(lincs$INCHIKEY,toxcast$INCHIKEY),"GSID"]

# Remove lincs compounds that dit not match any ToxCast compounds
lincs = subset(lincs,!is.na(TC_ID))

#Export to Rdata object
chemicals = lincs[,]
save(chemicals,file=file.path(DATA_DIR,"/intermediate_files/chemicals.Rdata"))

# save a copy to calculon as well
system(paste("scp ",file.path(DATA_DIR,"/intermediate_files/chemicals.Rdata")," calculon:",file.path(CALC_DATA_DIR,"boolGE"),sep=""))
