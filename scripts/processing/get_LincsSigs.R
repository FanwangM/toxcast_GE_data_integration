
#Path to the project directory - PLEASE CHANGE
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the data for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


# Libraries
library(plyr)
library(gdata)
library(jsonlite)

# Constants 
TIME= "24" #time at which the signatures were obtained after tratment (in hours)
DOSE = "10" #uM 

#Load chemicals, cellLines and gene annotations
load(file.path(CALC_DATA_DIR,"chemicals.Rdata"))
load(file.path(CALC_DATA_DIR,"boolGE/cellLines.Rdata"))
load(file.path(CALC_DATA_DIR,"boolGE/geneAnno_df.Rdata"))


# split compounds into batches of 100
batches = split(1:nrow(chemicals),1:10)


#initialize cellline counter
i_cell=0

# for each cellLine
l_ply(cellLines,function(cLine) {
    
    # increment cell counter
    i_cell<<-i_cell+1

    i_batch = 0
    # for each batch of compound
    lincsDF = llply(batches,function(b) {
          # extract the the top/down 50 compounds
          # for each signature corresponding to the
          # chemicals in the batch using custom url
          i_batch <<- i_batch + 1
          cat(cLine,": extraction for batch num",i_batch,"\n")
          url = URLencode(paste0('http://api.lincscloud.org/a2/siginfo?q={"pert_dose":"',DOSE,'","pert_itime":"',TIME,' h","pert_id":{"$in":["',paste0(chemicals$LINCS_ID[b],collapse='","'),'"]},"cell_id":"',cLine,'"}&f={"sig_id":1,"distil_cc_q75":1,"is_gold":1,"dn50_lm":1,"up50_lm":1}&l=1000&user_key=0cf1d0ca5ada162b28c379a0326ba2d3'))
          rd = try(readLines(url,warn=F),silent=T)
          df =  try(fromJSON(rd),silent=T)
          while(class(rd) == "try-error" | class(df) == "try-error") {
              cat("API's busy...\ntrying again in a few minutes\n")
              Sys.sleep(60)
              rd = try(readLines(url,warn=F),silent=T)
              df = try(fromJSON(rd),silent=T)
          }

          return(df)
      })

    # construct data frame from results of the batch
    lincsDF = do.call('rbind',lincsDF)
	
	
	# initialize chem counter
	i_chem = 0
	
    # for each compound
    sigs = llply(chemicals$LINCS_ID, function(chem) {
            
		#increment chem counter
		i_chem <<- i_chem+1
		#display message every 300 chemicals
		if(i_chem %% 300 == 0) 
                    cat(cLine,"(",i_cell,") - compound num:",i_chem,"\n")

        ## Fetch signatures corresponding to chemical
		chemData = lincsDF[grep(chem,lincsDF$sig_id),]
		
		# if result is not empty for this cell line
                if(nrow(chemData) > 0) {


                ## FIND BEST COMPOUND SIGNATURE
                # if multiple signatures for same compound and exp. conditions
                    if(nrow(chemData) > 1) {
                    # take reproducible and self-connected signatures if available
                        gold_sig = which(chemData$is_gold == T)
                        if(length(gold_sig) > 0)
                            chemData=chemData[gold_sig,]
                    # then take the one with the best spearman correlation across samples (to be confirmed)
                        chemData=chemData[which.max(chemData$distil_cc_q75),]
                    }

                    ## RETRIEVE GENE SYMBOLS FROM PROBESETS IDS
                    # for each up-regulated probeset_id
                    chemData$UP=list(llply(unlist(chemData$up50_lm),function(pbsetId) {
					     # retrieve gene symbol from gene annotation DF
						return(geneAnno.df[geneAnno.df$pr_id == pbsetId,"pr_gene_symbol"])
                    }))
                    # for each down-regulated probeset_id
                    chemData$DOWN=list(llply(unlist(chemData$dn50_lm),function(pbsetId){
						# same procedure
						return(geneAnno.df[geneAnno.df$pr_id == pbsetId,"pr_gene_symbol"])
                    }))
					
					# remove columns that are not needed anymore
                    chemData=chemData[,!(colnames(chemData) %in% c("_id","up50_lm","dn50_lm"))]
					
                    return(chemData)
                }
            })
		# save signature and associated genes in R object
        save(sigs,file=file.path(CALC_DATA_DIR,"boolGE",paste(cLine,"_sigs.Rdata",sep="")))
    })










 
