library(gficf)
library(plyr)
library(ggplot2)

plotbar = function(df.pred) {
  df.pred$lab = unlist(sapply(strsplit(x = rownames(df.pred),split = "_",fixed = T),function(x) x[1]))
  df.pred$type = data$embedded$type[match(df.pred$predicted.class,data$embedded$sample)]
  
  stat = ddply(df.pred,c("lab","predicted.class"),summarise,n=length(predicted.class))
  s = ddply(df.pred,c("lab"),summarise,n=length(predicted.class))
  stat$perc = stat$n/s$n[match(stat$lab,s$lab)]
  stat$predicted.class[stat$perc<0.05] = "Other"
  
  
  stat = ddply(stat,c("lab","predicted.class"),summarise,n=sum(n))
  s = ddply(df.pred,c("lab"),summarise,n=length(predicted.class))
  stat$perc = stat$n/s$n[match(stat$lab,s$lab)]
  
  print(ggplot(stat,aes(x=lab,y=perc,fill=predicted.class)) + geom_bar(stat = "identity",width = .25) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab(""))
  return(stat)
  
}

plotbar.MCF7 = function(df.pred) {
  df.pred$lab = unlist(sapply(strsplit(x = rownames(df.pred),split = "_",fixed = T),function(x) x[2]))
  df.pred$type = data$embedded$type[match(df.pred$predicted.class,data$embedded$sample)]
  
  stat = ddply(df.pred,c("lab","predicted.class"),summarise,n=length(predicted.class))
  s = ddply(df.pred,c("lab"),summarise,n=length(predicted.class))
  stat$perc = stat$n/s$n[match(stat$lab,s$lab)]
  stat$predicted.class[stat$perc<0.05] = "Other"
  
  
  stat = ddply(stat,c("lab","predicted.class"),summarise,n=sum(n))
  s = ddply(df.pred,c("lab"),summarise,n=length(predicted.class))
  stat$perc = stat$n/s$n[match(stat$lab,s$lab)]
  
  print(ggplot(stat,aes(x=lab,y=perc,fill=predicted.class)) + geom_bar(stat = "identity",width = .25) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab(""))
  
  return(stat)
  
}

plotPPV = function(df.pred) {
  require(precrec)
  require(pROC)
  df.pred$lab = unlist(sapply(strsplit(x = rownames(df.pred),split = "_",fixed = T),function(x) x[1]))
  df.pred$lab[df.pred$lab%in%"KPL1"] = "MCF7"
  df.pred$predicted.class[df.pred$predicted.class%in%"KPL1" & (df.pred$lab%in%"MCF7" | df.pred$lab%in%"KPL1")] = "MCF7"
  df.pred$tp = df.pred$predicted.class==df.pred$lab
  df.pred <- df.pred[order(df.pred$class.prob,df.pred$tp,decreasing = T),]
  df.pred$score = 10^6 - (1:length(df.pred$class.prob))
  
  #precrec_obj <- evalmod(scores = df.pred$score, labels = df.pred$tp*1)
  #autoplot(precrec_obj)
  #aucs <- evalmod(scores = df.pred$score, labels = df.pred$tp, mode = 'aucroc')
  #message("RAUC = ",round(aucs$uaucs$aucs,3))
  
  pROC_obj <- roc(df.pred$tp*1,df.pred$score,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=F, max.auc.polygon=F, grid=T,
                  print.auc=TRUE, show.thres=TRUE)
  
  #sens.ci <- ci.se(pROC_obj)
  #plot(sens.ci, type="shape", col="lightblue")
  return(pROC_obj)
}

plot.repositioning = function(data) {
  ggplot() + geom_point(data = data$embedded,aes(x=X,y=Y),shape=46,color="gray") + theme_bw() + geom_point(data = data$embedded.predicted,aes(x=X,y=Y),color="red",shape=46)
}

plotPPV.MCF7 = function(df.pred) {
  require(precrec)
  require(pROC)
  df.pred$lab = "MCF7"
  df.pred$lab[df.pred$lab%in%"KPL1"] = "MCF7"
  df.pred$predicted.class[df.pred$predicted.class%in%"KPL1" & (df.pred$lab%in%"MCF7" | df.pred$lab%in%"KPL1")] = "MCF7"
  df.pred$tp = df.pred$predicted.class==df.pred$lab
  df.pred <- df.pred[order(df.pred$class.prob,df.pred$tp,decreasing = T),]
  df.pred$score = 10^6 - (1:length(df.pred$class.prob))
  
  #precrec_obj <- evalmod(scores = df.pred$score, labels = df.pred$tp*1)
  #autoplot(precrec_obj)
  #aucs <- evalmod(scores = df.pred$score, labels = df.pred$tp, mode = 'aucroc')
  #message("RAUC = ",round(aucs$uaucs$aucs,3))
  
  pROC_obj <- roc(df.pred$tp*1,df.pred$score,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=F, max.auc.polygon=F, grid=T,
                  print.auc=TRUE, show.thres=TRUE)
  
  #sens.ci <- ci.se(pROC_obj)
  #plot(sens.ci, type="shape", col="lightblue")
  return(pROC_obj)
}


## MCF7 pred
data = loadGFICF("~/Downloads/BRCA.atlas.gficfv2.NMF.gficf")
M = readRDS("~/work/current/BRCA_AIRC_paper/paper_git/RData/MCF7_Broad.rds")
data = gficf::scMAP(data = data,x = M,nt = 16,seed = 0,normalize = T)
data = classify.cells(data = data,classes = data$embedded$sample,k = 101,knn_method = "manhattan",knn_weights_fun = NULL)
data$embedded.predicted$sample = unlist(sapply(strsplit(x = rownames(data$embedded.predicted),split = "_",fixed = T),function(x) x[2]))
s = plotbar.MCF7(data$embedded.predicted)
sum(s$n[s$predicted.class!="Other"])
pROC_obj <- plotPPV.MCF7(data$embedded.predicted)

## Cell line pan cancer
data = loadGFICF("~/Downloads/BRCA.atlas.gficfv2.NMF.gficf")
M = readRDS(file = "~/Downloads/agiv.revev.scRNA.RAW.rds")
df = data.frame(symb=unique(rownames(M)),stringsAsFactors = F)
df = gficf:::symbolToEns(df = df,col = "symb",organism = "human",verbose = T)
df = subset(df,!is.na(ens))
df = subset(df,!symb%in%df$symb[duplicated(df$symb)])
M = M[rownames(M)%in%df$symb,]
rownames(M) = df$ens[match(rownames(M),df$symb)]
colnames(M) = gsub(pattern = ".",replacement = "-",x = colnames(M),fixed = T)

meta = read.delim(file = "~/Downloads/Metadata.txt",stringsAsFactors = F)
rownames(meta) = meta$NAME
meta = subset(meta,Cancer_type%in%"Breast Cancer")
meta$Cell_line = gsub(pattern = "_BREAST",replacement = "",x = meta$Cell_line,fixed = T)
ccl <- c("CAMA1","ZR751","T47D","MDAMB361","KPL1","MCF7","BT474","MDAMB436","BT549","HCC38","HDQP1")
meta = subset(meta,Cell_line%in%ccl)
M = M[,meta$NAME]
M = M[rowSums(M)>0,]
colnames(M) = paste0(meta$Cell_line,"_",colnames(M))
#M = M[,colSums(M)>=5000]

data = gficf::scMAP(data = data,x = M,nt = 16,seed = 0,normalize = T)
data = classify.cells(data = data,classes = data$embedded$sample,k = 101,knn_method = "manhattan",knn_weights_fun = NULL)
pROC_obj <- plotPPV(data$embedded.predicted)

s = plotbar(data$embedded.predicted)

## Dropseq MDAMB468 & CAL51
M1 = readRDS("~/Downloads/MDAMB468_dropseq.rds")
data = loadGFICF("~/Downloads/BRCA.atlas.gficfv2.NMF.gficf")
data = gficf::scMAP(data = data,x = M1,nt = 16,seed = 0,normalize = T)
data = classify.cells(data = data,classes = data$embedded$sample,k = 101,knn_method = "manhattan",knn_weights_fun = NULL)
pROC_obj <- plotPPV(data$embedded.predicted)
plotbar(data$embedded.predicted)
plot.repositioning(data)

M2 = readRDS("~/Downloads/CAL51_dropseq.rds")
data = gficf::scMAP(data = data,x = M2,nt = 16,seed = 0,normalize = T)
data = classify.cells(data = data,classes = data$embedded$sample,k = 101,knn_method = "manhattan",knn_weights_fun = NULL)
pROC_obj <- plotPPV(data$embedded.predicted)
plotbar(data$embedded.predicted)
plot.repositioning(data)

