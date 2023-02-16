library(gficf)
library(ggplot2)
library(scRNAseq)
library(ggpubr)
library(reshape2)
library(ggrastr)
library(plyr)
library(SingleCellExperiment)
library(psupertime)
library(zellkonverter)


gmxPathways = function (gmx.file) 
{
  df = read.delim(gmx.file)
  pathways <- vector(mode = "list",length = ncol(df))
  names(pathways) = colnames(df)
  for (i in 1:ncol(df)) {pathways[[i]] <- df[,i]}
  pathways <- lapply(pathways, function(x) x[!x%in%""])
  pathways
}

plot.gene.sets = function(scores.fle,scores.gficf,fle.coord,p){
  fle.coord$pathw = scores.fle[,p]
  fle.coord$type = "wot"
  tmp = fle.coord
  tmp$pathw = scores.gficf[,p]
  tmp$type = "gficf.scGSEA"
  fle.coord = rbind(fle.coord,tmp);rm(tmp);gc()
  fle.coord$pathw[fle.coord$pathw<0]=0
  p1 = ggplot(data = subset(fle.coord,type%in%"wot"),aes(x=x,y=y,color=pathw)) + geom_point_rast(shape=46,size=0.01) + scale_color_distiller(palette = "Spectral") + theme_bw() + xlab("") + ylab("") + coord_fixed()
  p2 = ggplot(data = subset(fle.coord,type%in%"gficf.scGSEA"),aes(x=x,y=y,color=pathw)) + geom_point_rast(shape=46,size=0.01) + scale_color_distiller(palette = "Spectral") + theme_bw() + xlab("") + ylab("") + coord_fixed()
  
  ggarrange(p1,p2)
  
}

corr.scores = function(scores.fle,scores.gficf,p=c("Cell.cycle","Pluripotency","Neural.identity","Trophoblast","Epithelial.identity","MEF.identity")) {
  library(corrplot)
  library(RColorBrewer)
  r = cor(x = as.matrix(scores.fle[,p]),as.matrix(scores.gficf[,p],method = "spearman"))
  corrplot(r, order = 'FPC', type = 'lower', diag = T,col = rev(COL2('RdBu', 10)),addCoef.col = 'white',tl.pos = "l")
}

data = gficf::loadGFICF("~/Downloads/Schiebinger_data.gficf")
fle.coord = read.delim("~/Downloads/fle_coords.txt")
rownames(fle.coord) = fle.coord$id
fle.coord = fle.coord[rownames(data$embedded),]
scores.fle = read.csv("~/Downloads/gene_set_scores.csv",stringsAsFactors = F)
rownames(scores.fle) = scores.fle$id

corr.scores(scores.fle,scores.gficf = data$scgsea$x)

# rescale by gene sets
x = data$scgsea$x
x = t(x)
x = (x - rowMeans(x)) / apply(x, 1, sd)
x = t(x)

p=c("Cell.cycle","Pluripotency","Neural.identity","Trophoblast","Epithelial.identity","MEF.identity")
for (i in 1:length(p)){
  g = plot.gene.sets(scores.fle = scores.fle,scores.gficf = x,fle.coord = fle.coord,p = p[i])
  ggsave(filename = paste("~/Dropbox/GFICF_v2_manuscript/tmp/GFICF.wot.comparison",p[i],"pdf",sep = "."),plot = g,width = 6,height = 4)
}

