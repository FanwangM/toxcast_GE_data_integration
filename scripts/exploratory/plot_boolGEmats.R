

## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved 
## 3/11/14
## INPUT: boolean gene expression matrices
## OUTPUT: pdf with heatmaps


#Path to the project directory 
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'


# Libraries
library(plyr)
library(gplots)

# booleanGE matrices path and prefix
bGE_path = file.path(DATA_DIR,"boolGE")

# get list of booleanGE matrices files
files = system(paste("ls ",bGE_path,sep=""),intern=T)
files = files[grep("boolean_geneExp",files)]
l_ply(files,function(f) {
    cat(f,"\n")
    cline = unlist(strsplit(f,"boolean_geneExp_|.Rdata"))[2]
    load(file.path(bGE_path,f))
    df=get(ls()[grep("_df",ls())])
    if(ncol(df) > 0) {
        df = as.matrix(df)
        mode(df) = "numeric"
        png(file.path(PROJECT_DIR,"results/boolGE",paste(cline,"_htmap_clustered.png",sep="")))
        heatmap.2(df,trace="none",dendrogram="none",col=c("red","white","blue"))
        dev.off()
    }
})


