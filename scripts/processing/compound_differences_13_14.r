library(plyr)
library(gdata)

DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

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

toxcast = unique(toxcast)
dsstox = unique(dsstox)


dsstox=dsstox[-which(dsstox$INCHIKEY %in% toxcast$INCHIKEY),]

full_dat13 = read.xls("data/toxcast_2013/ToxCast_Generic_Chemicals_2013_12_10.xlsx",h=T)

# TOXCAST VS FULL DATASET13
which(is.na(match(toxcast$GSID,full_dat13$DSSTox_GSID))) 
# only 1 compound difference

# DSSToX vs DATASET13
dsstox$GSID2 = gsub("^","DSSTox_GSID_",dsstox$GSID)
which(is.na(match(dsstox$GSID2,full_dat13$DSSTox_GSID))) 
# again only 1 compound not matching original set
# which is very weird...




