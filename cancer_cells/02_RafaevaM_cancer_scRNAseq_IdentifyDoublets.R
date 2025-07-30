
#in this script, outputs from initial QC (01_RichinaV_CITEseq_HEALTHY_QC.R) are used. Doublets are identified


mydir <- "/scicore/home/bentires/GROUP/michal/rafaeva__maria/GSE223309_2023//human_data/starsolo_mapping_workflow"    # CHANGE!! working dir
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
  library(BiocSingular)
  library(tidyverse)
})

#====
# STEP 1: upload the fileterd data after QC and do lib normalization on the transcriptomic part
sce <- readRDS("Rdata/sce_demultiplexed_cleaned.rds")


library(scuttle)
#quick clustering for normalization
clust.trans <- quickCluster(sce) 
table(clust.trans)
#12clusters

sf.deconv <- calculateSumFactors(sce, cluster = clust.trans)
summary(sf.deconv)
#        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.1928  0.6610  0.8891  1.0000  1.2353  4.7687

sf.lib <- librarySizeFactors(sce); summary(sf.lib)

sizeFactors(sce) <- sf.deconv  #assign these factors to sce
sce <- logNormCounts(sce)

pdf("plots/libsize_trans.pdf", width = 4, height = 4)
plot(sf.lib, sf.deconv, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=clust.trans)
abline(a=0, b=1, col="red")
dev.off()


#PCA o tom HVGs
var.stats <- modelGeneVar(sce)
top.var <- getTopHVGs(sce, n=3000)
sce <- denoisePCA(sce, technical=var.stats, subset.row = top.var )
reducedDimNames(sce)
dim(reducedDim(sce,"PCA"))

#equip sce with the dimensioanl reduction
#sce <- runPCA(sce) #on log data!!
set.seed(101)
sce <- runTSNE(sce, dimred = "PCA")
set.seed(110)
sce <- runUMAP(sce, dimred = "PCA", ncomponents=5, n_neighbors=5)


#osca clustering


library(bluster)
set.seed(110001) #OVERCLUSTER
nn.clusters2 <- clusterCells(sce, use.dimred="PCA",  BLUSPARAM=SNNGraphParam(k=8, type="rank", cluster.fun="walktrap"))
table(nn.clusters2)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#738  587 1841 1448 1122 1046  268  381  553  647  843  330  324  484  151  150 
#17   18   19   20   21 
#65  114  278  533   24


sce$clusters.pca <- factor(nn.clusters2)

plotReducedDim(sce, "TSNE", colour_by="cell.line", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE)

pdf("plots/DimRedPlots_clusters_I.pdf", width = 20, height = 15)
par(mfrow=c(2,2))
gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="cell.line", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  plotReducedDim(sce, "PCA", colour_by="clusters.pca", text_colour = "red",text_by="clusters.pca",text_size = 3),
  plotReducedDim(sce, "TSNE", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  ncol = 2
)
dev.off()


#SOME CLUSTERS (SMALL ONES COULD BE DOUBLETS)
#Check doublets among the clusters
library(scDblFinder)
#https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
set.seed(132)
sce <- scDblFinder(sce,clusters = sce$clusters.pca)
table(sce$scDblFinder.class,sce$clusters.pca )
#         1    2    3    4    5    6    7    8    9   10   11   12   13   14
#singlet  650  548 1763 1365 1058  985    2  371  548  641  803  317  298  463
#doublet   88   39   78   83   64   61  266   10    5    6   40   13   26   21

#.         15   16   17   18   19   20   21
#singlet  146  149   62  106  260  528   18
#doublet    5    1    3    8   18    5    6

#sce2 <- scDblFinder(sce) #let them do the clustering
#table(sce2$scDblFinder.class,sce2$clusters.pca )

#dbl.out <- findDoubletClusters(sce,sce$clusters.pca) #needs to have clusters, this takes a little long
#it's checking similarity between each triplet of the clusters

#chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE, nmads = 3)]
#table(nn.clusters2)[chosen.doublet]



se <-  sumCountsAcrossCells(sce, sce$clusters.pca,  average=FALSE) #cretes SummarizedExp
assayNames(se) <- "counts"

keep <- rowSums(edgeR::cpm(assays(se)$counts) > 5) > 0 #maybe should use filterByExpr
sum(keep)
#14709 genes left
se <- se[keep,]
se <- logNormCounts(se)

GV.se <- modelGeneVar(se)
#Warning message:
# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
# collapsing to unique 'x' values
top.GV.se <- getTopHVGs(GV.se, prop=0.25)

z <- assays(se)$logcounts
ztop <- z[top.GV.se [1:100],] #100 most varying genes

library(pheatmap)
pdf("plots/heatmap_psudobulks_transcriptome_I.pdf", width = 8, height = 9)
pheatmap((ztop - rowMeans(ztop)), cex=.85)#
dev.off()


#save only cells annotated as singlets
table(sce$scDblFinder.class)
#singlet doublet 
#11081     846 

sce <- sce[,sce$scDblFinder.class ==  "singlet"]


#saveRDS(sce,"Rdata/sce_doublets_removed.rds")

