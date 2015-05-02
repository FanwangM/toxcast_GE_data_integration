
library(jsonlite)

CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

url="http://api.lincscloud.org/a1/geneinfo?q={%22pr_pool_id%22:{%22$regex%22:%22epsilon%22}}&f={%22pr_gene_symbol%22:1,%22pr_id%22:1}&l=1000&user_key=lincsdemo"
rd = try(readLines(url,warn=F),silent=T)
geneAnno.df = fromJSON(rd)
geneAnno.df = geneAnno.df[,-1]
save(geneAnno.df,file=file.path(CALC_DATA_DIR,"/boolGE/geneAnno_df.Rdata"))
