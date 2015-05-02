## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved
## 24/11/2014
## IN : Training Set + Selected Features + Test Set
## OUT : Likelihood of toxicity for the test set


### Initialisations ##########

#Path to the project directory
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


# Libraries
library(plyr)
library(randomForest)

# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
#TOX = "livrCarcino"
TOX='BW_decrs'

### Random Forest ####################

RF_Classif = function(trainingSet,testSet,selFeatures,nDS) {

	# Reduce training set to features selected by the MRMR
	# Generate levels according to the variable type
	# 0,1 for bioactivity and fingerprint
	# -1,0,1 for gene expression
	train = data.frame(lapply(1:length(selFeatures), function(i) {
						levels = unlist(ifelse((length(grep(paste(TOX,"|Tox_|FP_",sep=""),names(selFeatures)[i])) > 0),list(c(0,1)),list(c(-1,0,1))))
						return(factor(trainingSet[,selFeatures[i]],levels=levels))
						}))
	colnames(train) = names(selFeatures)

	# Same procedure for the test set
	test = data.frame(lapply(1:length(selFeatures), function(i) {
					levels = unlist(ifelse(length(grep(paste(TOX,"|Tox_|FP_",sep=""),names(selFeatures)[i])) > 0,list(c(0,1)),list(c(-1,0,1))))
					return(factor(testSet[,selFeatures[i]],levels=levels))
					}))
	colnames(test) = names(selFeatures)
	#tox_labels = test$liverCarcino
	test = subset(test,select=-c(BW_decrs))

	# Down-Sample negative class
	nontox.ind =  which(train[,TOX] == 0)
	tox.ind =  which(train[,TOX]== 1)

	pred = llply(1:nDS, function(i) {
           if(length(nontox.ind) > length(tox.ind)) {
		nontox.ind2 = sample(nontox.ind,length(tox.ind))
                tox.ind2 = tox.ind
            }
            else {
                tox.ind2 = sample(tox.ind,length(nontox.ind))
                nontox.ind2 = nontox.ind
            }

            train2 = train[c(nontox.ind2,tox.ind2),]
		#print(colnames(train2))

		# Training phase
		rf = randomForest(BW_decrs ~ ., data=train2, importance=T, ntree=300, na.action = na.omit)

		# Testing phase
		return(pred=list(predict(rf,test,type="vote",na.action=na.omit),imp=rf$importance))
	})
	return(pred)
}
