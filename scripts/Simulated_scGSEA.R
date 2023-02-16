require(gficf)
require(ggplot2)
require(reshape2)
require(corrplot)

# load data
sim.counts <- readRDS(file = "~/Downloads/Simulated_scRNAseq.rds")
sim.groups <- readRDS(file = "~/Downloads/simulated_groups.rds")
h.sets <- gficf:::gmtPathways(gmt.file = "~/Downloads/h.all.v7.4.symbols.gmt",convertToEns = F,convertHu2Mm = F)
sets.to.use <- 1:24
sets.and.groups <- data.frame(set=sets.to.use,group=paste("Group",rep(1:4,each=6),sep=""),set_name=names(h.sets)[sets.to.use],stringsAsFactors = F)

# Perform scGSEA
data = gficf::gficf(M = sim.counts,filterGenes = F,normalize = T)
data = gficf::runPCA(data = data,dim = 3)
data = gficf::runReduction(data = data,nt = 8,reduction = "umap",n_neighbors = 150)
data$embedded = cbind(data$embedded,sim.groups)
gficf::plotCells(data = data,colorBy = "Group") + xlab("UMAP1") + ylab("UMAP2")

h.sets = h.sets[names(h.sets)%in%sets.and.groups$set_name]
data = gficf::runScGSEA(data = data,geneID = "symbol",pathway.list = h.sets,nsim = 10000,nt = 31,verbose = T,nmf.k = 10,rescale = "byGS")#,use.for.nmf = "raw")

#### GS activity scores heatmap 
f = data$embedded[order(data$embedded$Group),]
df = melt(as.matrix(data$scgsea$x))
df$value[df$value< 0] = 0
df$Var1 = factor(as.character(df$Var1),levels = f$Cell)
df$Var2 = factor(as.character(df$Var2),levels = rev(sets.and.groups$set_name))
ggplot(data = df,aes(x=Var1,y=Var2,fill=value)) + geom_tile() + theme_bw() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + scale_fill_distiller(palette = "Blues",direction = 1)

### Overlapp between pathways
genes = unique(as.character(unlist(h.sets)))
M = Matrix::Matrix(data = 0,nrow = length(genes),ncol = length(h.sets))
rownames(M) = genes; colnames(M) = names(h.sets)
for (i in 1:ncol(M)) {M[h.sets[[i]],i] = 1}
jc = Matrix(data = 1- as.matrix(dist(t(M),method = "binary")),sparse = T)
colnames(jc) = rownames(jc) = paste0("GS",1:ncol(jc))

corrplot(as.matrix(jc),method = 'square', type = 'lower', diag = FALSE,is.corr = F,col = COL1("Reds"))

