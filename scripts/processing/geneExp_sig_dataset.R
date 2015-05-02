## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 13/10/14
## INPUT: Chemical Dataset
## OUTPUT: Gen expression DATaset with boolean variable for each signature_direction
## RUN ON CALCULON


### Initialisations ##########

#Path to the project directory 
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


# Libraries
library(plyr)
library(gdata)

# Load chemical dataset (name:chemicals,type:df)
load(file.path(CALC_DATA_DIR,"/chemicals.Rdata"))

# Load Cell-lines used in LINCS
load(file.path(CALC_DATA_DIR,"/boolGE/cellLines.Rdata"))

# get signature for each cell line
cLine_sigObjects  = system(paste("ls ",file.path(CALC_DATA_DIR,"boolGE")," | grep _sigs",sep=""),intern=T)


### Main loop ####################

# For each cell line
l_ply(cLine_sigObjects, function(cLine_obj) {

    # Load signatures (name:cLine_obj_[name of the cellLine],type:list of df)
    load(file.path(CALC_DATA_DIR,"boolGE",cLine_obj))
    cLine = unlist(strsplit(cLine_obj,"_"))[1]

    
    # Convert into data.frame
    sigs = do.call('rbind',sigs)
    
    # Create boolean gene expression data frame
    # with rows being the chemicals
    # and columns begin the gene postfixed by "up" or "down"
    boolGE_df=data.frame(rep(0,nrow(chemicals))) 
    rownames(boolGE_df) = chemicals$LINCS_ID
	
    # for each chemical in the chemical datasets
    i=0
    l_ply(rownames(boolGE_df), function(chem) {
        i <<- i+1
        if(i%%200 == 0)
            cat("iteration",i,"out of",nrow(boolGE_df),"- cellLine:",cLine,"\n")
        
        # get instances of the chemical in the gene expression dataset
        i_chem = grep(chem,sigs$sig_id)

        if(length(i_chem) > 0) {
            # for each gene signature matching the instances above
            # use the "regulation putter" function to put 1/-1 in the matrix
            up_sigs = sigs[i_chem,"UP"] 
            l_ply(unlist(up_sigs),function(geneSig) {
				#if gene signature does not exist in the boolean GE df
				if(!(geneSig %in% colnames(boolGE_df)))            
					boolGE_df[,geneSig]<<-rep(0,nrow(chemicals)) #then create it
				boolGE_df[chem,geneSig]<<-1
			})
            dn_sigs = sigs[i_chem,"DOWN"]
            l_ply(unlist(dn_sigs), function(geneSig) {
					if(!(geneSig %in% colnames(boolGE_df)))            
						boolGE_df[,geneSig]<<-rep(0,nrow(chemicals)) #then create it
					boolGE_df[chem,geneSig]<<-(-1)
			})
		}
	})


    # remove first column from matrix (used only to create the dataset)
    boolGE_df = boolGE_df[,-1]
    assign(paste(cLine,"_df",sep=""),boolGE_df)
    save(list=paste(cLine,"_df",sep=""),file=file.path(CALC_DATA_DIR,"boolGE",paste("boolean_geneExp_",cLine,".Rdata",sep="")))

})

# heatmaps
#library(gplots)
#heatmap.2(as.matrix(HA1E_df),trace="none",col=c("red","white","blue"))


