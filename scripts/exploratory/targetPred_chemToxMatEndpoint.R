
# load toxMat endpoints (name:toxmat, type:matrix)
load('~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/boolean_toxCast_toxEndpoints.Rdata')

# load chemicals (name:,type:df)
load('~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/chemicals.Rdata')

# get structures for chemicals in the tox mat endpoint
smi =  chemicals[match(rownames(toxmat),chemicals$LINCS_ID),'SMILES']

# set working dir to PIDGIN dir
setwd('~/Dropbox/ucc_az/PIDGIN/')

# write smiles in file in Pidgin dir
writeLines(smi,'tmp.smi')

# standardize according to PIDGIN protocol
system('~/ChemAxon/JChem/bin/standardize tmp.smi -c StandMoleProt.xml > tmp_st.smi')

# run target prediction tool (raw scores)
system('python predict_raw.py tmp_st.smi')

# run target prediction tool (thrshold predicitons)
system('python predict_binary_heat.py a tmp_st.smi')


# rename target prediction results
system('mv out_results.txt ~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/targetPred_rawScores_toxmat.txt')
system('mv out_results_binary_heat.txt ~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/targetPred_binaryResults_toxmat.txt')


# remove smiles files
system('rm tmp*.smi')

# open target prediciton files
targetScores = read.table('~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/targetPred_rawScores_toxmat.txt', h=F, sep='\t', comment.char='', quote='')
targetPred = read.table('~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/targetPred_binaryResults_toxmat.txt', h=F, sep='\t', comment.char='', quote='')


# rename columns according to chemicals in the analysis
colnames(targetScores) = c('Target','Uniprot',rownames(toxmat))
colnames(targetPred) = c('Target','Uniprot',rownames(toxmat))


# rewrite on file
write.table(targetScores,'~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/targetPred_rawScores_toxmat.txt',sep='\t',row.names=F,col.names=T,quote=F)
write.table(targetPred,'~/Dropbox/ucc_az/toxCast_lincs_integration/data/intermediate_files/targetPred_binaryResults_toxmat.txt',sep='\t',row.names=F,col.names=T,quote=F)

# melt binary prediction into two-column array
library(reshape2)
twoCols = melt(targetPred,ids=c('Name','Uniprot'))
twoCols = subset(twoCols,value==1,select=c(Name,variable)) #retain only predicted target/chemical pairs
#rename columns for clarity
colnames(twoCols) = c('Target','Compound')
write.table(twoCols,'two_columns_targetPred.txt',sep='\t',quote=F,row.names=F,col.names=T)

#END - communicate results to Marijke 


