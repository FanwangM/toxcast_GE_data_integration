## Supervisor : Dr Andreas Bender
## All rights reserved
## 26/11/2014
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
library(e1071)

# Cell line for gene expression dataset
CLINE = "MCF7"

#In-vivo Tox Outcome
#TOX = "livrCarcino"
TOX='BW_decrs'


### SVM ####################

SVM_Classif = function(trainingSet,testSet,selFeatures,nDS) {

	# Reduce training set to features selected by the MRMR
    train = trainingSet[,selFeatures]

	# Same procedure for the test set
	test = testSet[,selFeatures]
	test = subset(test,select=-c(BW_decrs))

	# For gene expression variables, split variables into up/down
	GE_vars = grep(paste(TOX,"|FP_|Tox21_",sep=""),colnames(train),invert=T)
	l_ply(GE_vars,function(i) {
	   gene = colnames(train)[i]
	   colnames(train)[i] <<- paste(gene,"_UP",sep="")
	   train[,paste(gene,"_DN",sep="")] <<- 0
	   train[which(train[,i] < 0),paste(gene,"_DN",sep="")] <<- 1
	   train[which(train[,i] < 0),i] <<- 0
	})
	l_ply(GE_vars,function(i) {
	   gene = colnames(test)[i]
	   colnames(test)[i] <<- paste(gene,"_UP",sep="")
	   test[,paste(gene,"_DN",sep="")] <<- 0
	   test[which(test[,i] < 0),paste(gene,"_DN",sep="")] <<- 1
	   test[which(test[,i] < 0),i] <<- 0

	})

	# Down-Sample larger class
	nontox.ind =  which(train[,TOX] == 0)
	tox.ind =  which(train[,TOX] == 1)

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
		svm = svm(BW_decrs~ ., data=train2, gamma=1, cost=1,probability=T, na.action = na.omit)

		# Testing phase
		return(predict(svm,test,probability=T, decision.values = T,na.action=na.omit))
	})
	return(pred)
}
