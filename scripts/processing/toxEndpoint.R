
library(reshape2)
library(RMySQL)
library(plyr)



CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

#chemicals
load(file.path(CALC_DATA_DIR,"/chemicals.Rdata"))

# load animal toxicity matrix 
load(file.path(CALC_DATA_DIR,'/datasets/bool_animalTox.Rdata'))
toxmat = table(animaltox$chid,as.character(animaltox$effect))

# remove effects with no matching chemical
toxmat = toxmat[, which(colSums(toxmat) > 0)]

# remove column with "nothing - effect"
toxmat = toxmat[,-1]

# replace count by hit
toxmat[toxmat > 0] = 1 

# replace toxCast id by Lincs ID
idx = match(rownames(toxmat),chemicals$TC_ID)
rownames(toxmat) = chemicals[idx,'LINCS_ID']
     
save(toxmat,file=file.path(CALC_DATA_DIR,'/datasets/boolean_toxCast_toxEndpoints.Rdata'))

