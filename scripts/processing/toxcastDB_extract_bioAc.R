
## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved
## 13/01/15
## EXTRACT IN-VITRO BIO.Ac DATA FROM THE TOXCAST DB (CALCULON)


library(reshape2)
library(RMySQL)
library(plyr)

####################
## Activity boolean
####################

CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

# chemicals (overlaping chems between Toxcast/Lincs)
load(file.path(CALC_DATA_DIR,"/chemicals.Rdata"))

drv = dbDriver("MySQL")
con = dbConnect(drv,user='az338',pass='pleasechange',dbname='toxcast')

# match chemicals GSID with sample IDs in the database
sampleIDs = llply(as.character(chemicals$TC_ID),function(gsid) {
    return(dbGetQuery(con,paste(
         "select sa_gsid, sa_sample_id from sample
          where sa_gsid = '",gsid,"'"
         ,sep='')))
})
df.sampleID = do.call('rbind',sampleIDs)

# Get hit calls (binary activities) and ac50s for each sample ID
hitC = llply(unique(df.sampleID$sa_sample_id),function(spid) {
    return(dbGetQuery(con,paste(
        "select spid ,l5.aeid, hitc, modl, modl_ga from level5 l5, level4 l4
         where l5.l4id = l4.l4id
         and l4.spid = '", spid, "'"
        ,sep='')))
})
df.hitC = do.call('rbind',hitC)

# Match sample ID by GSID
idx =  match(df.hitC$spid,df.sampleID$sa_sample_id)
df.hitC$GSID = df.sampleID[idx,'sa_gsid']

# Match GSID to Lincs id via corresponding toxcast ids
idx = match(df.hitC$GSID,chemicals$TC_ID)
df.hitC$LINCS_ID = chemicals[idx,'LINCS_ID']

# Extract assay name for each assay ids in results
aenm = llply(unique(df.hitC$aeid),function(id) {
             return(dbGetQuery(con,paste(
             #"select aeid, assay_component_endpoint_name
             "select * from assay_component_endpoint
             where aeid = '",as.character(id),"'"
             ,sep='')))
})
df.aenm = do.call('rbind',aenm)

# Match assay endpoint id with assay endpoint name
idx = match(df.hitC$aeid,df.aenm$aeid)
df.hitC$aenm = df.aenm[idx,'assay_component_endpoint_name']

# extract and save assay annotations
write.csv(df.aenm,file = file.path(CALC_DATA_DIR,'/datasets/assay_annotation.csv'),quote=F,row.names=F)


######################
##   HIT CALL MATRIX
######################
# Create boolean hit call matrix (chemicals vs assays)
# there should be at least 65 assays
# aggregation is set to max so that if compound is active
# at least once in a an assay, then it will be considered active

#hitcMat = acast(df.hitC,LINCS_ID ~ aenm,value.var = 'hitc',fun.aggregate=max)
hitcMat = acast(df.hitC, GSID  ~ aenm,value.var = 'hitc',fun.aggregate=max)
hitcMat[which(hitcMat == -Inf)] = NA

# Hit calls of -1 correspond to experiments
# with < 4 concentrations and do not correspond
# to active chemicals, therefore set to 0
hitcMat[which(hitcMat == -1)] = 0


# remove assays with no data point at all by counting number of
# values (not NA) in each assay
# updt : it seems that there is at least one data point in each assays
#valCounts = apply(hitcMat,2,function(x) length(which(!is.na(x))))
#sort(valCounts)
save(hitcMat,file = file.path(CALC_DATA_DIR,'/datasets/boolean_toxCast_hitCalls.Rdata'))



###################################################
## UNFILTERED AC50 MATRIX (project with Ain & Avid)
###################################################

ac50_mat = acast(df.hitC, LINCS_ID ~ aenm, value.var = 'modl_ga', fun.aggregate=median)
save(ac50_mat,file = file.path(CALC_DATA_DIR,'/datasets/ac50_toxCast_notFiltered.Rdata'))

###############################
## AC50 MATRIX (tox prediction)
###############################

## to be done :

# filter out hit calls of 0 since the estimated ac50 in
# this case does not correspond to a true active chemical
# for that particular assay


