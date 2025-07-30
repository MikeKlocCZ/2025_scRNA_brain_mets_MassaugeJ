

mydir <- "/scicore/home/bentires/GROUP/michal/rafaeva__maria/GSE223309_2023/human_data/starsolo_mapping_workflow"    # CHANGE!! working dir
setwd(mydir)
## Load packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(umap)
  library(DropletUtils)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(BiocParallel)
  library(BiocSingular)
  library(tidyverse)
  library(pheatmap)
})

#====
#ONLY CLEAN DATA NOW
# STEP 1: upload the fileterd data after QC and doublets removal, also upload the se with pseudobulks
sce <- readRDS("Rdata/sce_doublets_removed.rds")


#PCA o tom HVGs
var.stats <- modelGeneVar(sce)
top.var <- getTopHVGs(sce, n=3000)
sce <- denoisePCA(sce, technical=var.stats, subset.row = top.var )
reducedDimNames(sce)

#equip sce with the dimensioanl reduction
#sce <- runPCA(sce) #on log data!!
set.seed(101)
sce <- runTSNE(sce, dimred = "PCA")
set.seed(110)
sce <- runUMAP(sce, dimred = "PCA", ncomponents=5, n_neighbors=5)


#osca clustering


library(bluster)
set.seed(111) #
nn.clusters2 <- clusterCells(sce, use.dimred="PCA",    BLUSPARAM=SNNGraphParam(k=15, type="rank", cluster.fun="walktrap"))
table(nn.clusters2)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
#895 1330 2056 1695 1073  668  528  296  695 1276  258  112  122   77 

sce$clusters.pca <- nn.clusters2

pdf("plots/DimRedPlots_clusters_FINAL.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="clusters.pca", text_colour = "red",text_by="clusters.pca",text_size = 3),
  plotReducedDim(sce, "TSNE", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  ncol = 3
)
dev.off()

pdf("plots/DimRedPlots_clusters_FINAL_cellLine.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="cell.line" ),
  plotReducedDim(sce, "TSNE", colour_by="cell.line",add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="cell.line", add_legend = FALSE),
  ncol = 3
)
dev.off()



pdf("plots/DimRedPlots_clusters_FINAL_sampleType.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()



## BATCH CORRECTION BETWEEN EXPERIMENTS
library(batchelor)
set.seed(101)
f.out <- fastMNN(sce, batch=sce$Experiment, subset.row=top.var,  k = 25)
str(reducedDim(f.out, "corrected"))

# Transfer the correction results to the main spe object
reducedDim(sce, "fastMNN") <- reducedDim(f.out, "corrected")

set.seed(105)
sce <- runTSNE(sce, dimred="fastMNN")
set.seed(101)
sce <- runUMAP(sce, dimred="fastMNN", ncomponents=5, n_neighbors=10)


pdf("plots/FastMNN_batchCorrection.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()


# CLUSTER BASED ON BATCH CORRECTION


library(bluster)
set.seed(110) 
nn.clusters2 <- clusterCells(sce, use.dimred="fastMNN",  BLUSPARAM=SNNGraphParam(k=20, type="rank", cluster.fun="walktrap"))
table(nn.clusters2)
#     1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#801  122  469 1742  605  450  498 1112  501  267  489  182  477 1951  534  349 
#17   18   19   20   21 
#43   19  138  165  167 

sce$clusters.batchcor <- factor(nn.clusters2)

pdf("plots/DimRedPlots_clusters_FINAL_batchCorrected.pdf", width = 15, height = 8)
par(mfrow=c(1,2))
gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="clusters.batchcor", text_colour = "red", text_by="clusters.batchcor",text_size = 3,add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="clusters.batchcor", text_colour = "red", text_by="clusters.batchcor",text_size = 3,add_legend = FALSE),
  ncol = 2
)
dev.off()



#STEP 3: plotting the heatmaps

se <-  sumCountsAcrossCells(sce, sce$clusters.batchcor,  average=FALSE) #cretes SummarizedExp
assayNames(se) <- "counts"

keep <- rowSums(edgeR::cpm(assays(se)$counts) > 5) > 0 #maybe should use filterByExpr
sum(keep)
# 15476
se <- se[keep,]
se <- logNormCounts(se)
GV.se <- modelGeneVar(se)
#Warning message:
# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
# collapsing to unique 'x' values
top.GV.se <- getTopHVGs(GV.se, prop=0.25)

z <- assays(se)$logcounts
ztop <- z[top.GV.se [1:100],] #100 most varying genes

rownames(ztop) <- rowData(sce[rownames(ztop),])$SYMBOL

pdf("plots/pseudobulks_HM_transcriptome.pdf", width = 7, height = 9)
pheatmap((ztop - rowMeans(ztop)),fontsize=6)#
dev.off()



#saveRDS(sce,"Rdata/sce_ready_for_annotation.rds")
#saveRDS(se,"Rdata/se_ready_for_annotation.rds")
