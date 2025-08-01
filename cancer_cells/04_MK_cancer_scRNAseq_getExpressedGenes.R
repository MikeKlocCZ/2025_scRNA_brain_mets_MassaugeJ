
##
#  Extract expressed genes in the cell lines
#  from Masague paper, GSE223309

#Use Bioconductor 3.21 (R 4.5)
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(version = "3.21" ,force = TRUE)


pkg <- c("SingleCellExperiment",  "edgeR", "scuttle", "bluster",
         "scater", "scran", "umap", "tidyverse", "cowplot", "batchelor",
         "BiocParallel", "BiocSingular", "pheatmap")

pkg.install <- setdiff(pkg, rownames(installed.packages()))
if (length(pkg) > 0) {
  BiocManager::install(pkg.install)
}
## Load packages
suppressPackageStartupMessages({
  invisible(lapply(pkg, library, character.only = TRUE))
})


#====
#ONLY CLEAN DATA NOW
# STEP 1: upload the fileterd data after QC and doublets removal, also upload the se with pseudobulks
sce <- readRDS("Rdata/sce_ready_for_annotation.rds")

#MDA clusters 8, 14, 17, 9, 16, 20, 3, 10, 19
#HCC 13, 15, 12, 2, 21, 18, 1, 11, 7, 4, 5, 6


sce.MDA <- sce[,sce$clusters.batchcor %in% c("8", "14", "17", "9", "16", "20", "3", "10", "19")]
sce.HCC <- sce[,sce$clusters.batchcor %in% c("13", "15", "12", "2", "21", "18", "1", "11", "7", "4", "5", "6")]


#MDA
keep_feature <- rowSums(counts(sce.MDA) >  5) > 10
sum(keep_feature) #3263
sce.MDA <- sce.MDA[keep_feature,]

rowData(sce.MDA) |> names()
expressed.genes <-  rowData(sce.MDA)$SYMBOL[rowData(sce.MDA)$SYMBOL != ""]
 
write.csv(expressed.genes,"../TMEinteractions/MDAbr_expressed_genes.csv", row.names=FALSE)

#HCC
keep_feature <- rowSums(counts(sce.HCC) >  5) > 10
sum(keep_feature) #4418
sce.HCC <- sce.HCC[keep_feature,]

rowData(sce.HCC) |> names()
expressed.genes <-  rowData(sce.HCC)$SYMBOL[rowData(sce.HCC)$SYMBOL != ""]

write.csv(expressed.genes,"../TMEinteractions/HCCbr_expressed_genes.csv", row.names=FALSE)

