## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 24/11/2014
# IN : training set
# OUT : selected variables 



### Initialisations ##########

#Path to the project directory 
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


# Libraries
library(mRMRe)

# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
#TOX = "livrCarcino"
TOX='BW_decrs'


### Feature Selection ###########

MRMR_FtSel = function(trainDF, nF) {

	# change column classes to ordered factor
	# gene expressions were absolute-transformed to give same
	# importance to 1 and -1 

	rnames = rownames(trainDF)
	trainDF = data.frame(lapply(trainDF, function(x) ordered(abs(x))))
	rownames(trainDF) = rnames

	# identify tox column
	toxIndex = grep(TOX,colnames(trainDF))


	# Feature selection
	FS_data = mRMR.data(data = trainDF)
	FS = mRMR.classic(FS_data,target_indices=c(toxIndex),feature_count=nF)

	# Selected features + tox Column
	selFeatures = c(unlist(FS@filters),toxIndex)
	names(selFeatures) = c(colnames(trainDF)[unlist(FS@filters)],TOX)

	return(selFeatures)
	#save(selFeatures,file=file.path(CALC_DATA_DIR,"datasets/selFeatures.Rdata"))
}



