require(gficf)
require(psupertime)
require(Matrix)
require(ggplot2)
require(uwot)
require(SingleCellExperiment)

plot_all_genes_over_psupertime = function (psuper_obj,nrow,ncol,label_name = "Ordered labels",palette = "RdBu", plot_ratio = 1.25) 
{
  require(data.table)
  require(ggplot2)
  require(scales)
  make_col_vals <- function(y_labels, palette='RdBu')
  {
    require(RColorBrewer)
    require(grDevices)
    n_labels 	= length(levels(y_labels))
    max_col 	= 11
    if (n_labels==1) {
      col_vals 	= brewer.pal(3, palette)
      col_vals 	= col_vals[1]
    } else if (n_labels==2) {
      col_vals 	= brewer.pal(3, palette)
      col_vals 	= col_vals[-2]
    } else if (n_labels<=max_col) {
      col_vals 	= brewer.pal(n_labels, palette)
    } else {
      col_pal 	= brewer.pal(max_col, palette)
      col_vals 	= colorRampPalette(col_pal)(n_labels)
    }
    col_vals 	= rev(col_vals)
    
    return(col_vals)
  }
  
  proj_dt = psuper_obj$proj_dt
  beta_dt = psuper_obj$beta_dt
  x_data = psuper_obj$x_data
  ps_params = psuper_obj$ps_params
  top_genes = as.character(beta_dt$symbol)
  plot_wide = cbind(proj_dt, data.table(x_data[, top_genes,drop = FALSE]))
  plot_dt = melt.data.table(plot_wide, id = c("cell_id", "psuper","label_input", "label_psuper"), measure = top_genes, variable.name = "symbol")
  plot_dt[, symbol := factor(symbol, levels=top_genes)]
  col_vals = make_col_vals(plot_dt$label_input, palette)
  n_genes = length(top_genes)
  g = ggplot(plot_dt) + aes(x = psuper, y = value) + geom_point(size = 1,aes(colour = label_input)) + geom_smooth(se = FALSE, colour = "black") + scale_colour_manual(values = col_vals) +  scale_shape_manual(values = c(1, 16)) + scale_x_continuous(breaks = pretty_breaks()) + 
    scale_y_continuous(breaks = pretty_breaks()) + facet_wrap(~symbol, scales = "free_y", nrow = nrow, ncol = ncol) + theme_bw() + 
    theme(axis.text.x = element_blank(),legend.position = "top") + labs(x = "psupertime",y = "z-scored log2 expression", colour = label_name)
  return(g)
}

data = gficf::loadGFICF("~/Downloads/Aissa_dataset.gficf")
gficf::plotCells(data,colorBy = "time",pointSize = 1,pointShape = 19)

x = data$scgsea$x
x = x[,startsWith(x = colnames(x),prefix = "HALLMARK")]

y = sapply(strsplit(x = rownames(x),split = "-",fixed = T),function(x) x[1])
y = factor(as.character(y),levels = c("D0","D1","D2","D4","D9","D11"))
sce = SingleCellExperiment(list(logcounts=t(x)))
psuper_obj  = psupertime(sce, y, sel_genes='all',min_expression = 0)


g = plot_all_genes_over_psupertime(psuper_obj,nrow = 6,ncol = 8,label_name='Time')
print(g)
