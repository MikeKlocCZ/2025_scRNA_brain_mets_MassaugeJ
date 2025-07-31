
##
#  annotation using singleR and multiple reference data sets
# main reference is Allen Brain Atlas, however, since we are looking into metastatic brain
# we also use another reference like ImmGen to capture immune cells that might not be necessarily tissue resident
#
#  from Masague paper, GSE223309

#Use Bioconductor 3.21 (R 4.5)
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(version = "3.21" ,force = TRUE)


pkg <- c("SingleCellExperiment",  "edgeR", "scuttle","Seurat", 
         "scater", "scran", "umap", "tidyverse", "SingleR", "celldex",
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
# STEP 1: upload the fileterd data after QC and doublets removal, also upload the se with pseudobulks
sce <- readRDS("Rdata/sce_ready_for_annotation.rds")
rowData(sce)$ENSEMBL <- rownames(sce)
rownames(sce) <- rowData(sce)$SYMBOL


#This is the mouse brain atlas reference
ann.mouse.brain <- makeSummarizedExperimentFromLoom("input_data/mousebrain_reference/l5_all.agg.loom")

#explore and prepare for annotation
rownames(ann.mouse.brain) <- rowData(ann.mouse.brain)$Gene
names(assays(ann.mouse.brain))[1] <- "counts"
a <- table(rownames(ann.mouse.brain))
duplicated.names <- names(a[a > 1])
ann.mouse.brain <- ann.mouse.brain[!(rownames(ann.mouse.brain) %in% duplicated.names),]

logcounts <-  edgeR::cpm(ann.mouse.brain, log = TRUE)
assays(ann.mouse.brain,withDimnames=FALSE)$logcounts  <- logcounts 

ann.mouse.brain


#another annotation samples
ann.mouse.se <- MouseRNAseqData() 
ref.se.ImmGen <- ImmGenData()

#use celldex reference to annotate the cells
pred.each.cell.celldex <- SingleR(test = sce, ref = ann.mouse.se, assay.type.test=1, labels = ann.mouse.se$label.fine)
table(pred.each.cell.celldex$labels)

# Adipocytes                 aNSCs            Astrocytes  Astrocytes activated 
# 7                    44                    19                     3 
# B cells        Cardiomyocytes       Dendritic cells     Endothelial cells 
# 68                    19                     6                  2487 
# Ependymal          Erythrocytes           Fibroblasts Fibroblasts activated 
# 135                     5                   249                    87 
# Fibroblasts senescent          Granulocytes           Macrophages Macrophages activated 
# 4                   226                  1194                   112 
# Microglia   Microglia activated             Monocytes               Neurons 
# 11448                    76                   576                     7 
# NK cells                  NPCs      Oligodendrocytes                  OPCs 
# 70                   216                    91                     3 
# qNSCs               T cells 
# 276                    37 

tab.celldex <- table(Assigned=pred.each.cell.celldex$pruned.labels, Cluster=sce$clusters.batchcor)

pdf("FIGURES/04_Likelihood_cluster_cells_celldexMouseRNAseqData.pdf", width = 8, height = 8)
pheatmap(log2(tab.celldex+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()

#ImmGen
pred_ImmGen_fine <- SingleR(test = sce , ref = ref.se.ImmGen, 
                            assay.type.ref = "logcounts", labels = ref.se.ImmGen$label.main)

table(pred_ImmGen_fine$labels)
# B cells      B cells, pro         Basophils                DC Endothelial cells  Epithelial cells 
# 43                 3                 4               137              2536                 9 
# Fibroblasts               ILC       Macrophages         Microglia         Monocytes       Neutrophils 
# 682                47               797             12429               165               224 
# NK cells               NKT        Stem cells     Stromal cells           T cells               Tgd 
# 37                 8                 9               325                 9                 1 

tab_ImmGen_fine <- table(Assigned=pred_ImmGen_fine$pruned.labels, Cluster=sce$clusters.batchcor)

pdf("FIGURES/04_Likelihood_cluster_cells_ImmGen_main.pdf", width = 8, height = 6)
pheatmap(log2(tab_ImmGen_fine+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()

#mouse brain
pred.each.cell <- SingleR(test = sce, ref = ann.mouse.brain, assay.type.test=1, labels = ann.mouse.brain$ClusterName)
table(pred.each.cell$labels)

# ABC   ACBG   ACMB  ACNT1  ACNT2   ACOB  ACTE1  ACTE2  CBGRC   CHOR   COP1  ENMFB   ENT3   ENT4   ENT5 
# 119     17      7     12     23     11     31     43      2     16      3     44      1      1      1 
# ENT8  ENTG1  ENTG2  ENTG3  ENTG4  ENTG5  ENTG6  ENTG7   EPEN   EPMB   MGL1   MGL2   MGL3   MOL1   MOL2 
# 3     81     35     59     55      1     48     56    111      2   5114   3095   4795     59      2 
# MOL3  NFOL1 OBDOP2 OBNBL3    OEC    OPC   PER1   PER2   PER3   PVM1   PVM2   RGDG   RGSZ  SZNBL TEGLU6 
# 5      1      2     35    141     99     26      2     22    208    285      7     33     29      1 
# VECA   VECC   VECV  VLMC1  VLMC2  VSMCA 
# 5    263   2129      8     33    284 

tab <- table(Assigned=pred.each.cell$pruned.labels, Cluster=sce$clusters.batchcor)
pdf("FIGURES/04_Likelihood_cluster_cells_mouseBrain.pdf", width = 8, height = 8)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()


#****
# COMBINE THE ANNOT

#add cell names according to our annotation
clusters <- sce$clusters.batchcor

table(clusters)
clusters <- droplevels(clusters)

#CORRECTED ANNOTATION FROM CITEseq
# Classification each cell

levels(clusters)[levels(clusters) %in% c("2","9","8","14","27", "6", "19", "30", "26")] <- "Microglia"
levels(clusters)[levels(clusters) %in% c("3","1","31")] <- "Vascular endothelial cells, venous"
levels(clusters)[levels(clusters) %in% c("5","23")] <- "Perivascular macrophages"
levels(clusters)[levels(clusters) %in% c("18")] <- "Oligodendrocytes precursor cells"
levels(clusters)[levels(clusters) %in% c("28","33","35")] <- "Vascular smooth muscle cells, arterial"
levels(clusters)[levels(clusters) %in% c("4","24")] <-  "Olfactory ensheathing cells"
levels(clusters)[levels(clusters) %in% c("20")] <- "Vascular leptomeningeal cells"
levels(clusters)[levels(clusters) %in% c("42","16")] <- "Ependymal cells"
levels(clusters)[levels(clusters) %in% c("18")] <- "Neural stem cells"
levels(clusters)[levels(clusters) %in% c("32","22")] <- "Astrocytes"
levels(clusters)[levels(clusters) %in% c("21")] <- "Pericytes"
levels(clusters)[levels(clusters) %in% c("17","25")] <- "Oligodendrocytes"
levels(clusters)[levels(clusters) %in% c("34","7", "41")] <- "Vascular endothelial cells, capillary"
levels(clusters)[levels(clusters) %in% c("38")] <- "Subventricular zone\n radial glia-like cells"
levels(clusters)[levels(clusters) %in% c("39")] <- "Astrocyte-Bergmann glia"

  
#immuno
levels(clusters)[levels(clusters) %in% c("13",  "37")] <- "Neutrophils"
levels(clusters)[levels(clusters) %in% c("11")] <- "Monocytes"
levels(clusters)[levels(clusters) %in% c("10","40", "15")] <- "DCs"
levels(clusters)[levels(clusters) %in% c("20", "29")] <- "Fibroblasts"
levels(clusters)[levels(clusters) %in% c("36")] <- "B cells"
levels(clusters)[levels(clusters) %in% c("12")] <- "NK(T) cells"


sce$clusters <- clusters 

#save with annotation
#saveRDS(sce, "Rdata/sce_annotated_coarse.rds")

sce$location <- str_remove(sce$SampleName,"_.")
sce$cell.line.AND.condition <- paste(sce$cell.line,sce$SampleName, sep = "_") %>% str_remove("_\\d")
group.colors <- c(HCC1954_Stroma = "#333BFF", HCC1954_TME = "#9633FF", MDA231_Stroma  =  "#CC6600", MDA231_TME = "#E2FF33")


pdf("FIGURES/04_TSNE_annotation.pdf", width = 10, height = 10)
plotReducedDim(sce, "TSNE", colour_by="clusters", text_colour = "black", text_by="clusters",text_size = 4,add_legend = FALSE)
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
# 15047
se <- se[keep,]
se <- logNormCounts(se)
GV.se <- modelGeneVar(se)

#Warning message:
# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
# collapsing to unique 'x' values
top.GV.se <- getTopHVGs(GV.se, prop=0.025)

z <- assays(se)$logcounts
ztop <- z[top.GV.se,]


#out <- pheatmap((ztop - rowMeans(ztop)))#
#ordering.cols <- out$tree_col[["order"]]


pdf("FIGURES/04_pseudobulks_HM_transcriptome_clusterAnnot.pdf", width = 7, height = 14)
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
#6242  11223 

pdf("FIGURES/04_Barplot_Conditions_clusters_nonNormalized.pdf", width = 8, height = 5)
ggplot(data, aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

data.norm <- data
data.norm[data.norm$condition == "Stroma","value"] = 1/ max.cells["Stroma"]
data.norm[data.norm$condition == "TME","value"] = 1/ max.cells["TME"]

pdf("FIGURES/04_Barplot_Conditions_clusters.pdf", width = 8, height = 5)
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

pdf("FIGURES/04_Barplot_CellLines_clusters_nonNormalized.pdf", width = 8, height = 5)
ggplot(data, aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.off()

data.norm <- data
data.norm[data.norm$condition == "HCC1954","value"] = 1/ max.cells["HCC1954"]
data.norm[data.norm$condition == "MDA231","value"] = 1/ max.cells["MDA231"]

pdf("FIGURES/04_Barplot_CellLines_clusters.pdf", width = 8, height = 5)
ggplot(data.norm, aes(fill=condition, y=value, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions")+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.off()

