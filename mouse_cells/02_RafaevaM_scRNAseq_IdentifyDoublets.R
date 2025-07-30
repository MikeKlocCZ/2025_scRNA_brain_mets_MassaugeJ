
#in this script, outputs from initial QC (01_RichinaV_CITEseq_HEALTHY_QC.R) are used. Doublets are identified


mydir <- "/scicore/home/bentires/GROUP/michal/rafaeva__maria/GSE223309_2023/mouse_data/starsolo"    # CHANGE!! working dir
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
#26clusters

sf.deconv <- calculateSumFactors(sce, cluster = clust.trans)
summary(sf.deconv)
#        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.09937  0.47872  0.68306  1.00000  1.12298 15.01641

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
#1678  487   64  523 1561 2443  863 2395  360  749 1031  211  360  167  457  186 
#17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
#2246  149 1290  952  208   51   72  289  126  164  156   67   63   56   46   13 


sce$clusters.pca <- factor(nn.clusters2)

pdf("plots/DimRedPlots_clusters_I.pdf", width = 20, height = 8)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plotReducedDim(sce, "PCA", colour_by="clusters.pca", text_colour = "red",text_by="clusters.pca",text_size = 3),
  plotReducedDim(sce, "TSNE", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  plotReducedDim(sce, "UMAP", colour_by="clusters.pca", text_colour = "red", text_by="clusters.pca",text_size = 3,add_legend = FALSE),
  ncol = 3
)
dev.off()


#SOME CLUSTERS (SMALL ONES COULD BE DOUBLETS)
#Check doublets among the clusters
library(scDblFinder)
#https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
set.seed(132)
sce <- scDblFinder(sce,clusters = sce$clusters.pca)
table(sce$scDblFinder.class,sce$clusters.pca )
#             1    2    3    4    5    6    7    8    9   10   11   12   13   14
#singlet 1556  379   55  478 1551 2207  804 2013  356  744  942  187  160  162
#doublet  122  108    9   45   10  236   59  382    4    5   89   24  200    5

         #15   16   17   18   19   20   21   22   23   24   25   26   27   28
#singlet  455  186 1759  149 1219  765  208   49   70  266  126  162  155    6
#doublet    2    0  487    0   71  187    0    2    2   23    0    2    1   61

#         29   30   31   32
#singlet   61    1   20    0
#doublet    2   55   26   13

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
#14839 genes left
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
#ssinglet doublet 
#17251    2232

sce <- sce[,sce$scDblFinder.class ==  "singlet"]


#saveRDS(sce,"Rdata/sce_doublets_removed.rds")

