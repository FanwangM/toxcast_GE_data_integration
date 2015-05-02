## Author: Azedine Zoufir
## Supervisor : Dr Andreas Bender
## All rights reserved
## 29/11/14

library(plyr)
library(ggplot2)
library(reshape2)
library(gplots)

#Path to the project directory
PROJECT_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration/'
# Path to the directory containing the DATa for this project
DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

# Boolean matrix
load(file.path(DATA_DIR,'/intermediate_files/full_boolAct_mat.Rdata'))
boolAct_mat[which(is.na(boolAct_mat))] = -1
h_clust = heatmap.2(boolAct_mat)

boolAct_mat = as.data.frame(boolAct_mat)
boolAct_mat$Compound = rownames(boolAct_mat) 
boolDF = melt(boolAct_mat,id="Compound")
colnames(boolDF)[2:3] = c("Assay","Hit")

boolDF$Compound = factor(boolDF$Compound,levels=c(rownames(boolAct_mat)[h_clust$rowInd]))
boolDF$Assay = factor(boolDF$Assay,levels=c(colnames(boolAct_mat)[h_clust$colInd]))

ggplot(boolDF,aes(y=Compound,x=Assay,fill=as.factor(Hit))) +
geom_raster() +
scale_fill_manual(name="Hit", labels=c('NA','Inactive','Active'), values=c("white","beige","steelblue")) +
theme(axis.text.x=element_blank(),axis.text.y=element_blank())

# AC50
load(file.path(DATA_DIR,'/intermediate_files/full_glac50_mat.Rdata'))
glac50_mat[which(is.na(glac50_mat))] = -100
h_clust = heatmap.2(glac50_mat)
glac50_mat[which(glac50_mat == -100)] = NA

glac50_mat = as.data.frame(glac50_mat)
glac50_mat$Compound = rownames(glac50_mat) 
glac50DF = melt(glac50_mat,id="Compound")
colnames(glac50DF)[2:3] = c("Assay","log.ac50")

glac50DF$Compound = factor(glac50DF$Compound,levels=c(rownames(glac50_mat)[h_clust$rowInd]))
glac50DF$Assay = factor(glac50DF$Assay,levels=c(colnames(glac50_mat)[h_clust$colInd]))

ggplot(glac50DF,aes(y=Compound,x=Assay,fill=log.ac50)) +
  geom_raster() +
  scale_fill_gradient2(low="red",high="darkgoldenrod",na.value="white") +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())

#Animal tox
load(file.path(DATA_DIR,'/intermediate_files/bool_animalTox.Rdata'))
cttb = ddply(animaltox,.(chid,effect),summarize,n=length(effect))
cttb = acast(cttb,chid~effect,value.var="n")
cttb[which(is.na(cttb))] = 0
cttb[which(cttb > 0)] = 1
animaltox_mat= cttb[,-1]

h_clust = heatmap.2(animaltox_mat)

animaltox_mat = as.data.frame(animaltox_mat)
animaltox_mat$Compound = rownames(animaltox_mat) 
animaltox = melt(animaltox_mat,id="Compound")
colnames(animaltox)[2:3] = c("Tox_Endpoint","Hit")

animaltox$Compound = factor(animaltox$Compound,levels=c(rownames(animaltox_mat)[h_clust$rowInd]))
animaltox$Tox_Endpoint = factor(animaltox$Tox_Endpoint,levels=c(colnames(animaltox_mat)[h_clust$colInd]))

ggplot(animaltox,aes(y=Compound,x=Tox_Endpoint,fill=as.factor(Hit))) +
  geom_raster() +
  scale_fill_manual(name="Hit", labels=c('Inactive','Active'), values=c("beige","steelblue")) +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())

