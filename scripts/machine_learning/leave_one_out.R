## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 24/11/2014

# IN : full tox dataset 
# OUT : predictions for each compound
# PLEASE RUN ON CALCULON

### Initialisations ##########

# Libraries
library(plyr)

#Path to the project directory 
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_SCRIPTS_DIR='/home/az338/scripts/'

# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
#TOX = "livrCarcino"
TOX='BW_decrs'

### Leave one out ##########
# df = data
# nf = number of features 
# nDS = number of down sampling steps
# algo : machine learning algorithm to use ('rf' for randm forest, 'svm' for svm) 
# FP = logical indicating whether fingerprints should be included or not
# BioAc = logical indicating whether assay data (toxcast) should be included or not
# GE = logical indicating  whether gene expression data should be included or not  
lv1out = function(df,nF,nDS,algo,FP=T,BioAc=T,GE=T) {

	
	# if fingerprints not wanted in analysis 
	# remove variables from dataset
	if(!FP) df=df[,grep("FP_",colnames(df),invert=T)]
	
	# idem for bioactivity/assay variables
	if(!BioAc) df=df[,grep("Tox21_",colnames(df),invert=T)]
	
	# same for Gene Expression
	if(!GE) {
		tox = df[,TOX]
		df=df[,grep("FP_|Tox21_",colnames(df))]
		df[,TOX] = tox
	}
	# leave one out starts here
	out.pred = llply(1:nrow(df), function(i) {
	#out.pred = llply(1:5, function(i) {

		if((i %% 20) == 0) cat("Iter. ", i, "\n") 

		train = df[-i,]
		test = df[i,]
		

		# feature selection
		source(file.path(CALC_SCRIPTS_DIR,"Feature_Selection.R"))
		selcted_ft = MRMR_FtSel(train,nF)	
			

		# classification
		if(algo == 'rf') {
			source(file.path(CALC_SCRIPTS_DIR,"randomForest.R"))
			return(RF_Classif(train,test,selcted_ft,nDS))
		}
		else if(algo == 'svm') {
			source(file.path(CALC_SCRIPTS_DIR,"svm.R"))
			pred = try(SVM_Classif(train,test,selcted_ft,nDS))
			if(class(pred) == "try-error") pred = list(i,pred)
			return(pred)
		}
		else cat('Error : specified machine learning algorithm (',algo,') was not included in this script.')
	}) 
	return(out.pred)
}
