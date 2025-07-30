

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
sce <- readRDS("Rdata/sce_ready_for_annotation.rds")

pdf("plots/TSNE_afterBatcth.pdf", width = 10, height = 10)
plotReducedDim(sce, "TSNE", colour_by="cell.line", text_colour = "red", text_by="clusters.batchcor",text_size = 4, add_legend = TRUE)
plotReducedDim(sce, "TSNE", colour_by="detected", text_colour = "red", text_by="clusters.batchcor",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of detected genes")
plotReducedDim(sce, "TSNE", colour_by="sum", text_colour = "red", text_by="clusters.batchcor",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of total counts")
dev.off()

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
 
write.csv(expressed.genes,"TMEinteractions/MDAbr_expressed_genes.csv", row.names=FALSE)

#HCC
keep_feature <- rowSums(counts(sce.HCC) >  5) > 10
sum(keep_feature) #4418
sce.HCC <- sce.HCC[keep_feature,]

rowData(sce.HCC) |> names()
expressed.genes <-  rowData(sce.HCC)$SYMBOL[rowData(sce.HCC)$SYMBOL != ""]

write.csv(expressed.genes,"TMEinteractions/HCCbr_expressed_genes.csv", row.names=FALSE)


#MDA potential receptors
rownames(sce.MDA) <- rowData(sce.MDA)$SYMBOL

pdf("plots/TSNE_MDA_putative_receptors.pdf", width = 10, height = 10)
plotReducedDim(sce.MDA, "TSNE", colour_by="clusters.batchcor", text_colour = "red", text_by="clusters.batchcor",text_size = 4, add_legend = TRUE)
plotReducedDim(sce.MDA, "TSNE", colour_by="TGFBR1") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("TGFBR1 exprs")
plotExpression(sce.MDA, features="TGFBR1", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("TGFBR1 exprs")
plotReducedDim(sce.MDA, "TSNE", colour_by="TGFBR2") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("TGFBR2 exprs")
plotExpression(sce.MDA, features="TGFBR2", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("TGFBR1 exprs")
plotReducedDim(sce.MDA, "TSNE", colour_by="SDC2") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("SDC2 exprs")
plotExpression(sce.MDA, features="SDC2", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("SDC2 exprs")
plotReducedDim(sce.MDA, "TSNE", colour_by="RAMP1") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("RAMP1 exprs")
plotExpression(sce.MDA, features="RAMP1", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("RAMP1 exprs")
dev.off()

#MDA potential receptors
rownames(sce.HCC) <- rowData(sce.HCC)$SYMBOL

pdf("plots/TSNE_HCC_putative_receptors.pdf", width = 10, height = 10)
plotReducedDim(sce.HCC, "TSNE", colour_by="clusters.batchcor", text_colour = "red", text_by="clusters.batchcor",text_size = 4, add_legend = TRUE)
plotReducedDim(sce.HCC, "TSNE", colour_by="TGFBR1") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("TGFBR1 exprs")
plotExpression(sce.HCC, features="TGFBR1", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("TGFBR1 exprs")
plotReducedDim(sce.HCC, "TSNE", colour_by="SIRPA") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("SIRPA exprs")
plotExpression(sce.HCC, features="SIRPA", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("SIRPA exprs")

plotReducedDim(sce.HCC, "TSNE", colour_by="MMP13") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("MMP13 exprs")
plotExpression(sce.HCC, features="MMP13", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("MMP13 exprs")

plotReducedDim(sce.HCC, "TSNE", colour_by="IGF2R") +                                              # Change color brewer palette
  scale_color_distiller(palette = 2, direction = 1) + ggtitle("IGF2R exprs")
plotExpression(sce.HCC, features="IGF2R", 
               x="clusters.batchcor", colour_by="clusters.batchcor",add_legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 90))  + ggtitle("IGF2R exprs")
dev.off()


