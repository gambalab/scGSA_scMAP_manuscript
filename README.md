# GFICF version 2
# Single-cell gene set enrichment analysis and transfer learning for functional annotation of scRNA-seq data
#### M. Franchini<sup>1,3,+</sup>, S. Pellecchia<sup>1,+</sup>, G. Viscido<sup>1,+</sup>, G. Gambardella<sup>1,2,++</sup>

<sup>1</sup>Telethon Institute of Genetics and Medicine, Armenise/Harvard Laboratory of Integrative Genomics, Pozzuoli, Italy.  
<sup>2</sup>University of Naples Federico II, Department of Chemical, Materials and Industrial Engineering, Naples, Italy.  
<sup>3</sup>University of Naples Federico II, Department of Electrical Engineering and Information Technologies, 80125 Naples,Italy.  

## Abstract
Although an essential step, the functional annotation of cells often proves particularly challenging in the analysis of single-cell transcriptional data. Several methods have been developed to accomplish this task. However, in most cases, these rely on techniques initially developed for bulk RNA sequencing or simply make use of marker genes identified from cell clustering followed by supervised annotation. To overcome these limitations and automatise the process, we have developed two novel methods, the single-cell gene set enrichment analysis (scGSEA) and the single cell mapper (scMAP). scGSEA combines latent data representations and gene set enrichment scores to detect coordinated gene activity at single-cell resolution. scMAP uses transfer learning techniques to repurpose and contextualise new cells into a reference cell atlas. Using both simulated and real datasets, we show that scGSEA effectively recapitulates recurrent patterns of pathways’ activity shared by cells from different experimental conditions. At the same time, we show that scMAP can reliably map and contextualise new single cell profiles on a breast cancer atlas we recently released. Both tools are provided in an effective and straightforward workflow providing a framework to determine cell function and significantly improve annotation and interpretation of scRNA-seq data.

The full article [(Franchini, et al. 2022)](https://www.biorxiv.org/content/10.1101/2022.10.24.513476v1.full)

On-line tutorial can be found [HERE](https://github.com/gambalab/gficf)

Datasets used in the manuscript are available on figshare at following DOI: [10.6084/m9.figshare.22109414](https://figshare.com/articles/dataset/Single-cell_gene_set_enrichment_analysis_and_transfer_learning_for_functional_annotation_of_scRNA-seq_data/22109414)

## Figures

![alt text](https://github.com/gambalab/scGSA_scMAP_manuscript/blob/main/figures/Figure_01.png?raw=true)

<b>Figure 1 – Single cell gene set enrichment analysis overview and performances.</b> (A) GFICF package overview. (B) Single cell gene set enrichment analysis pipeline. (C) UMAP plot of 5,000 simulated cells grouped in four distinct groups. (D) Reconstructed activity of 24 simulated pathways across the 5,000 cells in (C). In the heatmap pathways are along rows while simulated cells along columns. Cells are ordered according to their group of origin. (E) Comparison between scGSEA pathway scores and signature scores originally computed by Schiebinger et al. on 251,203 single-cell profiles collected during differentiation stages. First row shows original gene set scores computed by Schiebinger et al. using wot phyton package. Second row shows gene set scores computed with scGSEA tool in the gficf R package. Each column represents a different gene set. Scores were plotted on the original FLE (force-directed layout) coordinates published by Schiebinger et al. (F) Spearman Correlation Coefficient (SCC) between scGSEA scores and wot package signature scores across the 251,203 single-cell transcriptional profiles in (E). (G) UMAP representation of 1,044 cells subject to eleven days of consecutive erlotinib treatment. Cells are color-coded according to sequenced day (i.e., 0, 1, 2, 4, 9, and 11 days). Single cell transcriptional profiles were normalized with gficf package. (H) EMT activity scores against inferred cell pseudo-time using the activity scores of 50 hallmark gene sets downloaded from MSigDB. Cells are color-coded as in (G). (I-K) Same as (H) but for wnt, cholesterol and fatty acid pathways respectively.
<hr/>

![alt text](https://github.com/gambalab/scGSA_scMAP_manuscript/blob/main/figures/Figure_02.png?raw=true)
<b>Figure 2 – Single cell mapping accuracy evaluation.</b> (A) UMAP representation of 35,276 cells from 32 breast cancer cell-lines using as cell subspace NMF. Cells are color-coded according to their cell-line of origin. Single cell transcriptional profiles were normalized with gficf package. (B) Same as (A) but using PC as cell sub-space (C) Accuracy evaluation of mapping method using cross-validation approach using as cell sub-space either NMF (orange) or PCA (green). Each boxplot display accuracy distribution in classifying the 32 cell lines but using as a training set the percentage of cells indicated on the x-axis. Accuracy is defined as the number of correctly classified cells over the total number of mapped cells. (D) Left plot; mapping of the MDAMB468 and CAL51 cells after they were re-sequenced with drop-SEQ technology. Right plot; accuracy of the mapping method in classifying re-sequenced MDAMB468 and CAL51 cells. NMF cell subspace is used. (E) Same as (D) but using PC sub-space. (F) Accuracy evaluation of the mapping method on 14,372 single cell transcriptomes sequenced with 10x Chromium method for MCF7 parental (P) cell line and three derived subclones (C1,C2,C3) using as cell sub-space either NMF (left) or PCA (right). (G) Performance of the mapping method in classifying 2,311 single cell transcriptomes sequenced using the 10x Chromium method from eleven distinct breast cell lines cell. Performances are reported in terms ROC curve and AUC is also displayed. (H) Same as (G) but using PC cell sub-space.
<hr/>

