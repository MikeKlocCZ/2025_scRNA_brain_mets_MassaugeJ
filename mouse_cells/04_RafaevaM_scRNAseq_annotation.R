## Set top directory
#Sys.setenv(DISPLAY="localhost:14.0")

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
  library(BiocParallel)
  library(BiocSingular)
  library(tidyverse)
})

#====
# STEP 1: upload the fileterd data after QC and doublets removal, also upload the se with pseudobulks
sce <- readRDS("Rdata/sce_ready_for_annotation.rds")
rowData(sce)$ENSEMBL <- rownames(sce)
rownames(sce) <- rowData(sce)$SYMBOL

library(loomR)
library(Seurat)
library(SeuratDisk)

#SINGLE CELL
#loom.mousebrain <- connect(filename = "mousebrain_reference/l5_all.loom", mode = "r")
#seurat.mousebrain <- as.Seurat(loom.mousebrain )
#sce.mousebrain <- as.SingleCellExperiment(seurat.mousebrain )
#saveRDS(sce.mousebrain , file = "mousebrain_reference/sce_mousebrain.rds")

#https://github.com/pachterlab/kb_python/issues/26
#loom.mousebrain.agg <- connect(filename = "mousebrain_reference/l5_all.agg.loom", mode = "r", skip.validate = TRUE)
ann.mouse.brain <- makeSummarizedExperimentFromLoom("mousebrain_reference/l5_all.agg.loom")

#explore and prepare for annotation
rownames(ann.mouse.brain) <- rowData(ann.mouse.brain)$Gene
names(assays(ann.mouse.brain))[1] <- "counts"
a <- table(rownames(ann.mouse.brain))
duplicated.names <- names(a[a > 1])
ann.mouse.brain <- ann.mouse.brain[!(rownames(ann.mouse.brain) %in% duplicated.names),]

logcounts <-  edgeR::cpm(ann.mouse.brain, log = TRUE)
assays(ann.mouse.brain,withDimnames=FALSE)$logcounts  <- logcounts 

ann.mouse.brain


library(celldex)
library(SingleR)
library(pheatmap)

#another annotation samples
ann.mouse.se <- MouseRNAseqData() 
ref.se.ImmGen <- ImmGenData()

pred.each.cell.celldex <- SingleR(test = sce, ref = ann.mouse.se, assay.type.test=1, labels = ann.mouse.se$label.fine)
table(pred.each.cell.celldex$labels)

tab.celldex <- table(Assigned=pred.each.cell.celldex$pruned.labels, Cluster=sce$clusters.batchcor)

pdf("plots/Likelihood_cluster_cells_celldexMouseRNAseqData.pdf", width = 8, height = 8)
pheatmap(log2(tab.celldex+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()

#ImmGen
pred_ImmGen_fine <- SingleR(test = sce , ref = ref.se.ImmGen, 
                            assay.type.ref = "logcounts", labels = ref.se.ImmGen$label.main)

table(pred_ImmGen_fine$labels)
tab_ImmGen_fine <- table(Assigned=pred_ImmGen_fine$pruned.labels, Cluster=sce$clusters.batchcor)
pdf("plots/Likelihood_cluster_cells_ImmGen_main.pdf", width = 8, height = 6)
pheatmap(log2(tab_ImmGen_fine+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()

#mouse brain
pred.each.cell <- SingleR(test = sce, ref = ann.mouse.brain, assay.type.test=1, labels = ann.mouse.brain$ClusterName)
table(pred.each.cell$labels)

#AABC   ACBG   ACMB  ACNT1  ACNT2   ACOB  ACTE1  ACTE2  CBGRC   CHOR   COP1 
#121     19      8     11     23     11     31     38      2     17      3 
#ENMFB   ENT3   ENT4   ENT5   ENT8  ENTG1  ENTG2  ENTG3  ENTG4  ENTG5  ENTG6 
#44      1      1      2      3     81     35     58     56      1     48 
#ENTG7   EPEN   EPMB  HYPEN   MGL1   MGL2   MGL3   MOL1   MOL2   MOL3  NFOL1 
#57    113      2      1   5094   3054   4656     51      2      5      1 
#OBDOP2 OBNBL3    OEC    OPC   PER1   PER2   PER3   PVM1   PVM2   RGDG   RGSZ 
#2     36    142     99     25      2     21    204    282      8     35 
#SZNBL TEGLU6   VECA   VECC   VECV  VLMC1  VLMC2  VSMCA 
#30      1      5    264   2147      7     31    260

tab <- table(Assigned=pred.each.cell$pruned.labels, Cluster=sce$clusters.batchcor)
pdf("plots/Likelihood_cluster_cells_mouseBrain.pdf", width = 8, height = 8)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()


#****
# COMBINE THE ANNOT

#add cell names according to our annotation
#remove 39
sce <- sce[,!(sce$clusters.batchcor == "39")]
clusters <- sce$clusters.batchcor


clusters <- droplevels(clusters)

#CORRECTED ANNOTATION FROM CITEseq
# Classification each cell

levels(clusters)[levels(clusters) %in% c("7","14","5","3","11", "24", "32")] <- "Microglia"
levels(clusters)[levels(clusters) %in% c("4","2","29")] <- "Vascular endothelial cells, venous"
levels(clusters)[levels(clusters) %in% c("10","18")] <- "Perivascular macrophages"
levels(clusters)[levels(clusters) %in% c("13")] <- "Oligodendrocytes precursor cells"
levels(clusters)[levels(clusters) %in% c("36","26","23")] <- "Vascular smooth muscle cells, arterial"
levels(clusters)[levels(clusters) %in% c("6","30")] <-  "Olfactory ensheathing cells"
levels(clusters)[levels(clusters) %in% c("22")] <- "Vascular leptomeningeal cells"
levels(clusters)[levels(clusters) %in% c("37")] <- "Ependymal cells"
levels(clusters)[levels(clusters) %in% c("25")] <- "Neural stem cells"
levels(clusters)[levels(clusters) %in% c("31")] <- "Astrocytes"
levels(clusters)[levels(clusters) %in% c("27")] <- "Pericytes"
levels(clusters)[levels(clusters) %in% c("19","28")] <- "Oligodendrocytes"
levels(clusters)[levels(clusters) %in% c("34")] <- "Neuroblasts"
levels(clusters)[levels(clusters) %in% c("1")] <- "Vascular endothelial cells, capillary"
levels(clusters)[levels(clusters) %in% c("38")] <- "Enteric glia"
levels(clusters)[levels(clusters) %in% c("16")] <- "Subventricular zone\n radial glia-like cells"
  
  
#immuno
levels(clusters)[levels(clusters) %in% c("12","33")] <- "Neutrophils"
levels(clusters)[levels(clusters) %in% c("8","15")] <- "Macrophages/Monocytes"
levels(clusters)[levels(clusters) %in% c("9","21")] <- "DCs"
levels(clusters)[levels(clusters) %in% c("20","16")] <- "Fibroblasts"
levels(clusters)[levels(clusters) %in% c("17")] <- "B cells"
levels(clusters)[levels(clusters) %in% c("35")] <- "NK(T) cells"

#levels(clusters)[levels(clusters) %in% c("39")] <- "remove"


sce$clusters <- clusters 


#saveRDS(sce, "Rdata/sce_annotated_coarse.rds")
sce$location <- str_remove(sce$SampleName,"_.")
sce$cell.line.AND.condition <- paste(sce$cell.line,sce$SampleName, sep = "_") %>% str_remove("_\\d")
group.colors <- c(HCC1954_Stroma = "#333BFF", HCC1954_TME = "#9633FF", MDA231_Stroma  =  "#CC6600", MDA231_TME = "#E2FF33")


pdf("plots/TSNE_annotation.pdf", width = 10, height = 10)
plotReducedDim(sce, "TSNE", colour_by="clusters", text_colour = "red", text_by="clusters",text_size = 4,add_legend = FALSE)
plotReducedDim(sce, "TSNE", colour_by="SampleName", text_colour = "black", text_by="clusters",text_size = 5,add_legend = TRUE)
plotReducedDim(sce, "TSNE", colour_by="cell.line", text_colour = "black", text_by="clusters",text_size = 5,add_legend = TRUE) + ggtitle("TME and stroma combined")
plotReducedDim(sce[,sce$cell.line == "MDA231"], "TSNE", colour_by="location", text_colour = "black", text_by="clusters",text_size = 5,add_legend = TRUE) + ggtitle("TME and stroma, only MDA231")
plotReducedDim(sce[,sce$cell.line == "HCC1954"], "TSNE", colour_by="location", text_colour = "black", text_by="clusters",text_size = 5,add_legend = TRUE) + ggtitle("TME and stroma, only HCC1945")
plotReducedDim(sce, "TSNE", colour_by="cell.line.AND.condition", text_colour = "black", text_by="clusters",text_size = 5,add_legend = TRUE) +  scale_color_manual(values=group.colors)
plotReducedDim(sce, "TSNE", colour_by="Experiment", text_colour = "black", text_by="clusters",text_size = 5,add_legend = TRUE)
dev.off()



se <-  sumCountsAcrossCells(sce, sce$clusters , average=FALSE) #cretes SummarizedExp
assayNames(se) <- "counts"

keep <- rowSums(edgeR::cpm(assays(se)$counts) > 5) > 0 #maybe should use filterByExpr
sum(keep)
# 15120
se <- se[keep,]
se <- logNormCounts(se)
GV.se <- modelGeneVar(se)

#Warning message:
# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
# collapsing to unique 'x' values
top.GV.se <- getTopHVGs(GV.se, prop=0.025)

z <- assays(se)$logcounts
#ztop <- z[top.GV.se [1:100],] #100 most varying genes
ztop <- z[top.GV.se,]


#out <- pheatmap((ztop - rowMeans(ztop)))#
#ordering.cols <- out$tree_col[["order"]]


pdf("plots/pseudobulks_HM_transcriptome_clusterAnnot.pdf", width = 7, height = 14)
pheatmap((ztop - rowMeans(ztop)),fontsize=5)#
dev.off()




# Compostion of the clusters
#for barplots we need a columns with ones
cluster <- sce$clusters
condition <- sce$SampleName |> str_remove("_\\d")
cell.line <- sce$cell.line
value <- rep(1,dim(sce)[2])
data <- data.frame(cluster,condition,value, cell.line)
# Stacked

table(cluster,condition)
max.cells <- table(condition)
# Stroma    TME 
#   6194  11057

pdf("plots/Barplot_Conditions_clusters_nonNormalized.pdf", width = 8, height = 5)
ggplot(data, aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

data.norm <- data
data.norm[data.norm$condition == "Stroma","value"] = 1/ max.cells["Stroma"]
data.norm[data.norm$condition == "TME","value"] = 1/ max.cells["TME"]

pdf("plots/Barplot_Conditions_clusters.pdf", width = 8, height = 5)
ggplot(data.norm, aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, both cell lines")

ggplot(data.norm[data.norm$cell.line == "MDA231",] , aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, MDA231")

ggplot(data.norm[data.norm$cell.line == "HCC1954",] , aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, HCC1954")


dev.off()


#for cell lines
# Compostion of the clusters
#for barplots we need a columns with ones
cluster <- sce$clusters
condition <- sce$cell.line
value <- rep(1,dim(sce)[2])
data <- data.frame(cluster,condition,value)
# Stacked

table(cluster,condition)
max.cells <- table(condition)
# HCC1954  MDA231 
# 9597    7568 

pdf("plots/Barplot_CellLines_clusters_nonNormalized.pdf", width = 8, height = 5)
ggplot(data, aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.off()

data.norm <- data
data.norm[data.norm$condition == "HCC1954","value"] = 1/ max.cells["HCC1954"]
data.norm[data.norm$condition == "MDA231","value"] = 1/ max.cells["MDA231"]

pdf("plots/Barplot_CellLines_clusters.pdf", width = 8, height = 5)
ggplot(data.norm, aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions")+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.off()

