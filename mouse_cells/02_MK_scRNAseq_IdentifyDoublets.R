

##
#  Dublet detection, mouse brain single-cell data 
#  from Masague paper, GSE223309

#Use Bioconductor 3.21 (R 4.5)
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(version = "3.21" ,force = TRUE)


pkg <- c("SingleCellExperiment",  "edgeR", "scuttle", "bluster",
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
# STEP 1: upload the fileterd data after QC and do lib normalization on the transcriptomic part
sce <- readRDS("Rdata/sce_demultiplexed_cleaned.rds")


#quick clustering for normalization by deconvolution
set.seed(121)
clust.trans <- quickCluster(sce) 
table(clust.trans)
#26clusters

sf.deconv <- calculateSumFactors(sce, cluster = clust.trans)
summary(sf.deconv)
#       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0880  0.4872  0.6920  1.0000  1.1212 14.4065 #

sf.lib <- librarySizeFactors(sce); summary(sf.lib)

sizeFactors(sce) <- sf.deconv  #assign these factors to sce
sce <- logNormCounts(sce)

pdf("FIGURES/02_libsize_trans.pdf", width = 4, height = 4)
plot(sf.lib, sf.deconv, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=clust.trans)
abline(a=0, b=1, col="red")
dev.off()


#PCA o tom HVGs
var.stats <- modelGeneVar(sce)
top.var <- getTopHVGs(sce, n=3000)
sce <- denoisePCA(sce, technical=var.stats, subset.row = top.var )
reducedDimNames(sce) #PCA added
dim(reducedDim(sce,"PCA")) # 19483     5

#equip sce with the dimensioanl reduction
set.seed(101)
sce <- runTSNE(sce, dimred = "PCA")
set.seed(110)
sce <- runUMAP(sce, dimred = "PCA", ncomponents=5, n_neighbors=5)

# Searching for doublets
#osca clustering
set.seed(110001) #OVERCLUSTER
nn.clusters2 <- clusterCells(sce, use.dimred="PCA",  BLUSPARAM=SNNGraphParam(k=8, type="rank", cluster.fun="walktrap"))
table(nn.clusters2)
#      1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21 
# 531  570  380  239 1003  829  191 1413  242 1000  508 2098 1753  208  382 1229  172 1558  853  864 1070 
# 22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38 
# 125  648  121   50  181   67  511  135   32   97  144   64   59   56   42   45   13 


sce$clusters.pca <- factor(nn.clusters2)

pdf("FIGURES/02_DimRedPlots_clusters_I.pdf", width = 20, height = 8)
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
#https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
set.seed(132)
sce <- scDblFinder(sce,clusters = sce$clusters.pca)
table(sce$scDblFinder.class,sce$clusters.pca )
#         1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
# singlet  485  518  179  213  996  783  191 1200  226  976  424 1680 1585  206  346 1187  167 1305  790
# doublet   46   52  201   26    7   46    0  213   16   24   84  418  168    2   36   42    5  253   63
# 
# 20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38
# singlet  857 1043  121  593  121   47  180   62  445  133   20   94  141   63   15    1   41   31    0
# doublet    7   27    4   55    0    3    1    5   66    2   12    3    3    1   44   55    1   14   13


se <-  sumCountsAcrossCells(sce, sce$clusters.pca,  average=FALSE) #cretes SummarizedExp
assayNames(se) <- "counts"

keep <- rowSums(edgeR::cpm(assays(se)$counts) > 5) > 0 #filter low expressed genes
sum(keep)
#14935 genes left
se <- se[keep,]
se <- logNormCounts(se)

GV.se <- modelGeneVar(se)
#Warning message:
# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
# collapsing to unique 'x' values
top.GV.se <- getTopHVGs(GV.se, prop=0.25)

z <- assays(se)$logcounts
ztop <- z[top.GV.se [1:100],] #100 most varying genes

#visual check on the pseudobulk level
pdf("FIGURES/02_heatmap_psudobulks_transcriptome_I.pdf", width = 8, height = 9)
pheatmap((ztop - rowMeans(ztop)), cex=.85)#
dev.off()


#save only cells annotated as singlets
table(sce$scDblFinder.class)
#singlet doublet 
#17465    2018 

sce <- sce[,sce$scDblFinder.class ==  "singlet"]
#saveRDS(sce,"Rdata/sce_doublets_removed.rds")

