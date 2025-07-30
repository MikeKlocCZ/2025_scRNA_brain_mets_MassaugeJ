

##
#  Proper clustering of a clean data, mouse brain single-cell data 
#. We'll also correct for Batch Effect using fast MNN
#  from Masague paper, GSE223309

#Use Bioconductor 3.21 (R 4.5)
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(version = "3.21" ,force = TRUE)


pkg <- c("SingleCellExperiment",  "edgeR", "scuttle", "bluster", "batchelor",
         "scater", "scran", "umap", "tidyverse", "cowplot", "scDblFinder",
         "BiocParallel", "BiocSingular", "pheatmap")

pkg <- setdiff(pkg, rownames(installed.packages()))
if (length(pkg) > 0) {
  BiocManager::install(pkg)
}
## Load packages
suppressPackageStartupMessages({
  invisible(lapply(pkg, library, character.only = TRUE))
})


#====
#ONLY CLEAN DATA NOW
# STEP 1: upload the fileterd data after QC and doublets removal, also upload the se with pseudobulks
sce <- readRDS("Rdata/sce_doublets_removed.rds")


#PCA o tom HVGs
var.stats <- modelGeneVar(sce)
top.var <- getTopHVGs(sce, n=3000)
sce <- denoisePCA(sce, technical=var.stats, subset.row = top.var )

#equip sce with the dimensional reduction on clean data
set.seed(101)
sce <- runTSNE(sce, dimred = "PCA")
set.seed(110)
sce <- runUMAP(sce, dimred = "PCA", ncomponents=5, n_neighbors=5)


#osca clustering
set.seed(111) #
nn.clusters2 <- clusterCells(sce, use.dimred="PCA",  BLUSPARAM=SNNGraphParam(k=15, type="rank", cluster.fun="walktrap"))
table(nn.clusters2)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21 
# 772 1444  943 2634  187  801  252 2041 1715  271  711 1733  177  204  160  423 1244  209  321  243  131 
# 22   23   24 
# 585   71  193 

sce$clusters.pca <- nn.clusters2

pdf("FIGURES/03_DimRedPlots_clusters_FINAL.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="clusters.pca", text_colour = "red",text_by="clusters.pca",text_size = 3),
  plotReducedDim(sce, "TSNE", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  ncol = 3
)
dev.off()

pdf("FIGURES/03_DimRedPlots_clusters_FINAL_cellLine.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="cell.line" ),
  plotReducedDim(sce, "TSNE", colour_by="cell.line",add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="cell.line", add_legend = FALSE),
  ncol = 3
)
dev.off()



pdf("FIGURES/03_DimRedPlots_clusters_FINAL_sampleType.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()


#Just Exp1

sce1 <- sce[,sce$Experiment == 1]
pdf("FIGURES/03_DimRedPlots_clusters_FINAL_sampleType_exp1.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce1, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce1, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce1, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()

sce2 <- sce[,sce$Experiment == 2]
pdf("FIGURES/03_DimRedPlots_clusters_FINAL_sampleType_exp2.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce2, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce2, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce2, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()


## BATCH CORRECTION BETWEEN EXPERIMENTS
set.seed(101)
f.out <- fastMNN(sce, batch=sce$Experiment, subset.row=top.var,  k = 25)
str(reducedDim(f.out, "corrected"))

# Transfer the correction results to the main spe object
reducedDim(sce, "fastMNN") <- reducedDim(f.out, "corrected")

set.seed(105)
sce <- runTSNE(sce, dimred="fastMNN")
set.seed(101)
sce <- runUMAP(sce, dimred="fastMNN", ncomponents=5, n_neighbors=10)


pdf("FIGURES/03_FastMNN_batchCorrection.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()


sce1 <- sce[,sce$Experiment == 1]
pdf("FIGURES/03_DimRedPlots_clusters_FINAL_sampleType_exp1_batchCor.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce1, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce1, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce1, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()

sce2 <- sce[,sce$Experiment == 2]
pdf("FIGURES/03_DimRedPlots_clusters_FINAL_sampleType_exp2_batchCor.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce2, "PCA", colour_by="SampleName" ),
  plotReducedDim(sce2, "TSNE", colour_by="SampleName",add_legend = FALSE),
  plotReducedDim(sce2, "UMAP", colour_by="SampleName", add_legend = FALSE),
  ncol = 3
)
dev.off()


#########
# CLUSTER BASED ON BATCH-CORRECTED DATA
set.seed(110) 
nn.clusters2 <- clusterCells(sce, use.dimred="fastMNN",  BLUSPARAM=SNNGraphParam(k=20, type="rank", cluster.fun="walktrap"))
table(nn.clusters2)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21 
# 297 3712 1530   78  545  560  300 3340 2086   60  128  106  197 1176   65   65   36  155   45   55   22 
# 22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42 
# 78  236  113   41  372  955  169   43   98  123   72   66  142   29   38   20   25  119   38   20  110 

sce$clusters.batchcor <- factor(nn.clusters2)

pdf("FIGURES/03_DimRedPlots_clusters_FINAL_batchCorrected.pdf", width = 15, height = 8)
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
# 15499
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

pdf("FIGURES/03_pseudobulks_HM_transcriptome.pdf", width = 7, height = 9)
pheatmap((ztop - rowMeans(ztop)),fontsize=6)#
dev.off()



#saveRDS(sce,"Rdata/sce_ready_for_annotation.rds")
#saveRDS(se,"Rdata/se_ready_for_annotation.rds")
