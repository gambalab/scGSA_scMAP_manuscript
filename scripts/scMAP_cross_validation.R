library(gficf)
library(ggplot2)
library(ggbeeswarm)
library(plyr)

# Download count matrix from https://figshare.com/ndownloader/files/28893384
M.raw = readRDS(file = "~/Downloads/RAW.UMI.counts.BC.cell.lines.rds")
M.raw = M.raw[Matrix::rowSums(M.raw)>0,]
subspace = "NMF"

stat = NULL
set.seed(0)

for (perc in c(10,20,30,40,50,60,70,80,90)) {
  print(perc)
  cells.test = NULL
  samples = unlist(sapply(strsplit(x = colnames(M.raw),split = "_",fixed = T),function(x) x[1]))
  u = unique(samples)
  for(i in u)
  {
    ix = samples %in% i
    n = sum(ix) - round(sum(ix)/100*perc)
    c = colnames(M.raw)[ix]
    cells.test = c(cells.test,sample(x = c,size = n))
  }
  
  M.test = M.raw[,cells.test]
  M.test = M.test[Matrix::rowSums(M.test)>0,]
  M.train = M.raw[,!(colnames(M.raw)%in%cells.test)]
  M.train = M.train[Matrix::rowSums(M.train)>0,]
  
  ccl.type = read.delim("~/ccl.classification.txt",stringsAsFactors = F)
  data = gficf::gficf(M = M.train,verbose = T,filterGenes = T,normalize = T);
  if (subspace == "PCA") {
    data = gficf::runPCA(data = data,dim = 50)
  } else {
    data = gficf::runNMF(data = data,dim = 100)    
  }
  data = gficf::runReduction(data = data,"umap",n_neighbors=15,nt = 16)
  data$embedded$sample = unlist(sapply(strsplit(x = colnames(data$gficf),split = "_",fixed = T),function(x) x[1]))
  data$embedded$cell.id = rownames(data$embedded)
  
  gficf::plotCells(data = data,colorBy = "sample",pointShape = 46,pointSize = .01)
  
  data = gficf::embedNewCells(data = data,x = M.test,nt = 4,verbose = T)
  data = classify.cells(data = data,classes = data$embedded$sample,k = 101,knn_method = "manhattan",knn_weights_fun = NULL)
  
  df = data$embedded.predicted
  df$lab = sapply(strsplit(x = rownames(df),split = "_",fixed = T),function(x) x[[1]])
  df$tp = df$lab == df$pred
  print(sum(df$tp) / nrow(df) *100)
  
  stat.gficf = ddply(df,"lab",summarise,ppv=sum(tp)/length(tp))
  stat.gficf$type = perc
  stat = rbind(stat,stat.gficf)
  gc()
}

s = ddply(stat,"type",summarise,mu=mean(ppv),med=median(ppv),sig=sd(ppv))

ggplot(data = stat,aes(x=as.character(type),y=ppv)) + geom_boxplot(outlier.colour =  NA,width=.25,fill="orange") + ylim(c(.95,1)) + theme_bw() 
