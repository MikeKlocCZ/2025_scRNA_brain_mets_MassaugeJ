
##
# Here we further look into the annotated cell types
# We perform Cell cycle analysis using cyclone from scran
# Identify marker genes for each cell type (internal annotation check)
# We also look into major cell clusters
#  from Masague paper, GSE223309

#Use Bioconductor 3.21 (R 4.5)
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(version = "3.21" ,force = TRUE)


pkg <- c("SingleCellExperiment",  "edgeR", "scuttle","Seurat",  "velociraptor",
         "scater", "scran", "umap", "tidyverse", "SingleR", "celldex", "bluster",
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
sce <- readRDS("Rdata/sce_annotated_coarse.rds")
colLabels(sce) <- sce$clusters

write.table(table(colLabels(sce)),"tables/cellTypeProps.txt")

# plot technical and cell cycle differences in the clusters
cyclin.genes <- grep("^Ccn[abde][0-9]$", rowData(sce)$SYMBOL)
cyclin.genes <- rownames(sce)[cyclin.genes]
cyclin.genes

pdf("FIGURES/05_CyclinsExprs.pdf", width = 10, height = 5)
plotHeatmap(sce, order_columns_by="label", 
            cluster_rows=FALSE, features=sort(cyclin.genes))
dev.off()


#cyclone
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)

pdf("FIGURES/05_cellcycle.pdf", width = 7, height = 7)
plot(assignments$score$G1, assignments$score$G2M,
     xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()

table(assignments$phases, colLabels(sce))
sce$cyclone.phase <- assignments$phases


#equip the annotated sce with the cyclone score
saveRDS(sce,"Rdata/sce_annotated_coarse.rds")

pdf("FIGURES/05_TSNE_technical_parameters.pdf", width = 10, height = 10)
plotReducedDim(sce, "TSNE", colour_by="detected", text_colour = "black", text_by="clusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of detected genes")
plotReducedDim(sce, "TSNE", colour_by="sum", text_colour = "red", text_by="clusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of total counts")
plotReducedDim(sce, "TSNE", colour_by="cyclone.phase", text_colour = "red", text_by="clusters",text_size = 4, add_legend = TRUE)
dev.off()

#############
#marker genes of the clusters
marker.info <- scoreMarkers(sce, colLabels(sce)) #length equal to tne number of clusters
cluster.names <- levels(sce$clusters)


pdf("FIGURES/05_MarkerGenesExprs.pdf", width = 10, height = 9)
for (ii in seq_along(cluster.names)){
  chosen <- marker.info[[ cluster.names[ii] ]]
  ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
  fig <- plotExpression(sce, features=head(rownames(ordered), n = 10), 
                 x="label", colour_by="label",add_legend = FALSE) +
                scale_x_discrete(guide = guide_axis(angle = 90)) + ggtitle(paste(cluster.names[ii],"top 10 marker genes"))
  print(fig)
}
dev.off()


#subclustering of the big clusters
# 1 microglia
# 2 Vascular endothelial cells 
# 3 Perivascular macrophages

table(sce$clusters)
sce$FinerClusters <- as.character(sce$clusters)


#MICROGLIA SUBSET
sce.Microglia <- sce[,colLabels(sce) %in% c("Microglia")]
set.seed(101) 
nn.clusters.Microglia <- clusterCells(sce.Microglia, use.dimred="fastMNN",  BLUSPARAM=SNNGraphParam(k=25, type="rank", cluster.fun="walktrap"))
table(nn.clusters.Microglia)
#  1    2    3    4    5    6    7    8    9   10   11   12   13 
#547 1703  868 1082 2158  524  435 1989  114 1185 1451  256   32 

#cluster Microgila.1 is the one cycling
sce.Microglia$subclusters <- paste0("Microglia.",nn.clusters.Microglia)

saveRDS(sce.Microglia,"Rdata/sceMicroglia_clusters.rds")
#sce.Microglia <- readRDS("Rdata/sceMicroglia_clusters.rds")

#ADD FINER CLUSTERS
sce$FinerClusters[colnames(sce) %in% colnames(sce.Microglia)]  <- sce.Microglia$subclusters
#plotReducedDim(sce, "TSNE", colour_by="FinerClusters", text_colour = "black", text_by="FinerClusters",text_size = 4,add_legend = TRUE)

###########################
#trajectories
dec.Microglia <- modelGeneVar(sce.Microglia)
top.hvgs.Microglia <- getTopHVGs(dec.Microglia, n=2000)
velo.out.Microglia <- scvelo(sce.Microglia, subset.row=top.hvgs.Microglia, assay.X="spliced")

sce.Microglia$velocity_pseudotime <- velo.out.Microglia$velocity_pseudotime
plotReducedDim(sce.Microglia,dim = "fastMNN", colour_by="velocity_pseudotime")

embedded.TSNE <- embedVelocity(reducedDim(sce.Microglia, "TSNE"), velo.out.Microglia)
grid.df <- gridVectors(sce.Microglia, embedded.TSNE, use.dimred = "TSNE")

embedded.UMAP <- embedVelocity(reducedDim(sce.Microglia, "UMAP"), velo.out.Microglia)
grid.df.UMAP <- gridVectors(sce.Microglia, embedded.UMAP, use.dimred = "UMAP")
#############################
#saveRDS(sce.Microglia,"Rdata/sceMicroglia_clusters.rds")
#sce.Microgila <- "Rdata/sceMicroglia_clusters.rds"

#interface to scVelo
#https://statbiomed.github.io/SingleCell-Workshop-2021/RNA-velocity.html

#reticulate::use_condaenv("sgcell", required = TRUE)
#library(reticulate)
#py_config()
#conda_list()
#use_condaenv("sgcell")
#plt <- import("matplotlib.pyplot", as = "plt")
#scv <- import ("scvelo")


pdf("FIGURES/05_TSNE_microglia_annotation.pdf", width = 10, height = 10)
plotReducedDim(sce.Microglia, "TSNE", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)
plotReducedDim(sce.Microglia, "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("both cell lines")
plotReducedDim(sce.Microglia[,sce.Microglia$cell.line == "MDA231"], "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("MDA231")
plotReducedDim(sce.Microglia[,sce.Microglia$cell.line == "HCC1954"], "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("HCC1954")

plotReducedDim(sce.Microglia, "UMAP", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) + xlim(-4.3,1.5) + ylim(-4.8,5)
plotReducedDim(sce.Microglia, "TSNE", colour_by="detected", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of detected genes")
plotReducedDim(sce.Microglia, "TSNE", colour_by="sum", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of total counts")
plotReducedDim(sce.Microglia, "TSNE", colour_by="Experiment", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Which experiment")
plotReducedDim(sce.Microglia[,sce.Microglia$Experiment == "1"], "TSNE", colour_by="cyclone.phase", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Cycle stage in Experiment 1")
plotReducedDim(sce.Microglia[,sce.Microglia$Experiment == "2"], "TSNE", colour_by="cyclone.phase", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Cycle stage in Experiment 2")
plotReducedDim(sce.Microglia, "TSNE", colour_by="cell.line", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Cycle stage in Experiment 2")
plotTSNE(sce.Microglia, colour_by="velocity_pseudotime") +
  geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, 
                                         xend=end.1, yend=end.2, colour=NULL), arrow=arrow(length=unit(0.05, "inches"))) 
plotReducedDim(sce.Microglia,  "UMAP", colour_by="velocity_pseudotime") +
  geom_segment(data=grid.df.UMAP, mapping=aes(x=start.1, y=start.2, 
                                         xend=end.1, yend=end.2, colour=NULL), arrow=arrow(length=unit(0.05, "inches")))   + xlim(-4.3,1.5) + ylim(-4.8,5)

dev.off()



#####
# bar plots
#
#

cluster <- sce.Microglia$subclusters
condition <- sce.Microglia$cell.line
location <- str_remove(sce.Microglia$SampleName,"_.")
Experiment <- sce.Microglia$Experiment
data <- data.frame(cluster,condition,Experiment,location)
# Stacked

table(cluster,condition)
max.cells <- table(condition)
# HCC1954  MDA231 
#

max.cells.Exp <- table(Experiment)
# 1     2 
#10431  1913 

max.cells.location <- table(location)
  #Stroma    TME 
#2454   9890 


data.norm <- data
data.norm$scaling.condition <-  data.norm$scaling.Experiment  <-  data.norm$scaling.location <- 1 

data.norm[data.norm$condition == "HCC1954",]$scaling.condition = 1/ max.cells["HCC1954"]
data.norm[data.norm$condition == "MDA231",]$scaling.condition = 1/ max.cells["MDA231"]

data.norm[data.norm$Experiment == 1,]$scaling.Experiment = 1/ max.cells.Exp[1]
data.norm[data.norm$Experiment == 2,]$scaling.Experiment = 1/ max.cells.Exp[2]


pdf("FIGURES/05_MG_Barplot_clusters.pdf", width = 8, height = 5)
ggplot(data.norm, aes(fill=location, y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, both cell lines")

ggplot(data.norm[data.norm$condition == "MDA231",], aes(fill=location , y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, Only MDA231") + scale_fill_manual(values=c("darkgreen",  "darkblue"))

ggplot(data.norm[data.norm$condition == "HCC1954",], aes(fill=location , y=scaling.condition, x=cluster)) +  
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, Only HCC1954") + scale_fill_manual(values=c("darkgreen",  "darkblue"))

ggplot(data.norm, aes(fill=condition, y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, cluster specificity wrt. cell lines")+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


ggplot(data.norm, aes(fill=Experiment, y=scaling.Experiment, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, differences in Experiments")+ scale_fill_manual(values=c("lightgreen",  "#56B4E9"))

dev.off()


#markers microglia
#marker genes of the clusters
colLabels(sce.Microglia) <- as.factor(sce.Microglia$subclusters)
marker.info.glia <- scoreMarkers(sce.Microglia, colLabels(sce.Microglia)) #length equal to tne number of clusters
cluster.names <- levels(colLabels(sce.Microglia) )




pdf("FIGURES/05_MG_MarkerGenesExprs.pdf", width = 10, height = 10)
for (ii in seq_along(cluster.names)){
  chosen <- marker.info.glia[[ cluster.names[ii] ]]
  ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
  fig <- plotExpression(sce.Microglia, features=head(rownames(ordered), n = 10), 
                        x="label", colour_by="label",add_legend = FALSE) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + ggtitle(paste(cluster.names[ii],"top 10 marker genes"))
  print(fig)
}
dev.off()




#VASCULAR ENDOTHELIAL
sce.VECs <- sce[,colLabels(sce) %in% c("Vascular endothelial cells, venous","Vascular endothelial cells, capillary")]
set.seed(101) 
nn.clusters.VECs <- clusterCells(sce.VECs, use.dimred="fastMNN",  BLUSPARAM=SNNGraphParam(k=25, type="rank", cluster.fun="walktrap"))
table(nn.clusters.VECs)
#   1    2    3    4    5    6    7    8 
#275 1102  117  145  172  304  181  116 
sce.VECs$subclusters <- paste0("Vascular endothelial cells.",nn.clusters.VECs)



#plotReducedDim(sce.VECs, "TSNE", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)

#ADD FINER CLUSTERS
sce$FinerClusters[colnames(sce) %in% colnames(sce.VECs)]  <- sce.VECs$subclusters




pdf("FIGURES/05_TSNE_VECs_annotation.pdf", width = 8, height = 7)
plotReducedDim(sce.VECs, "TSNE", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)
plotReducedDim(sce.VECs, "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("both cell lines")
plotReducedDim(sce.VECs[,sce.VECs$cell.line == "MDA231"], "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("MDA231")
plotReducedDim(sce.VECs[,sce.VECs$cell.line == "HCC1954"], "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("HCC1954")
plotReducedDim(sce.VECs, "TSNE", colour_by="detected", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of detected genes")
plotReducedDim(sce.VECs, "TSNE", colour_by="sum", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of total counts")
plotReducedDim(sce.VECs, "TSNE", colour_by="Experiment", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Which experiment")
plotReducedDim(sce.VECs[,sce.VECs$Experiment == "1"], "TSNE", colour_by="cyclone.phase", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Cycle stage in Experiment 1")
plotReducedDim(sce.VECs[,sce.VECs$Experiment == "2"], "TSNE", colour_by="cyclone.phase", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Cycle stage in Experiment 2")
plotReducedDim(sce.VECs, "TSNE", colour_by="cell.line", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("cell line")
dev.off()



#####
# bar plots VECs
#
#

cluster <- sce.VECs$subclusters
condition <- sce.VECs$cell.line
location <- str_remove(sce.VECs$SampleName,"_.")
Experiment <- sce.VECs$Experiment
data <- data.frame(cluster,condition,Experiment,location)
# Stacked

table(cluster,condition)
max.cells <- table(condition)
# HCC1954  MDA231 
#  1088    1324 

max.cells.Exp <- table(Experiment)
# 1     2 
#1030 1382



data.norm <- data
data.norm$scaling.condition <-  data.norm$scaling.Experiment  <-  data.norm$scaling.location <- 1 

data.norm[data.norm$condition == "HCC1954",]$scaling.condition = 1/ max.cells["HCC1954"]
data.norm[data.norm$condition == "MDA231",]$scaling.condition = 1/ max.cells["MDA231"]

data.norm[data.norm$Experiment == 1,]$scaling.Experiment = 1/ max.cells.Exp[1]
data.norm[data.norm$Experiment == 2,]$scaling.Experiment = 1/ max.cells.Exp[2]


pdf("FIGURES/05_VECs_Barplot_clusters.pdf", width = 7, height = 5)
ggplot(data.norm, aes(fill=location, y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, both cell lines")
ggplot(data.norm[data.norm$condition == "MDA231",], aes(fill=location , y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, Only MDA231") + scale_fill_manual(values=c("darkgreen",  "darkblue"))
ggplot(data.norm[data.norm$condition == "HCC1954",], aes(fill=location , y=scaling.condition, x=cluster)) +  
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, Only HCC1954") + scale_fill_manual(values=c("darkgreen",  "darkblue"))
ggplot(data.norm, aes(fill=condition, y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, cluster specificity wrt. cell lines")+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
ggplot(data.norm, aes(fill=Experiment, y=scaling.Experiment, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, differences in Experiments")+ scale_fill_manual(values=c("lightgreen",  "#56B4E9"))
dev.off()




#PERIVASCULAR MACRO
sce.PeriMacro <- sce[,colLabels(sce) %in% c("Perivascular macrophages")]
set.seed(101) 
nn.clusters.PeriMacro <- clusterCells(sce.PeriMacro, use.dimred="fastMNN",  BLUSPARAM=SNNGraphParam(k=30, type="rank", cluster.fun="walktrap"))
table(nn.clusters.PeriMacro)
#    1   2   3 
#188 374 219 
sce.PeriMacro$subclusters <- paste0("Perivascular macrophages.",nn.clusters.PeriMacro)


#ADD FINER CLUSTERS
sce$FinerClusters[colnames(sce) %in% colnames(sce.PeriMacro)]  <- sce.PeriMacro$subclusters




pdf("FIGURES/05_PeriMacro_VECs_annotation.pdf", width = 7, height = 6)
plotReducedDim(sce.PeriMacro, "TSNE", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)
plotReducedDim(sce.PeriMacro, "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("both cell lines")
plotReducedDim(sce.PeriMacro[,sce.PeriMacro$cell.line == "MDA231"], "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("MDA231")
plotReducedDim(sce.PeriMacro[,sce.PeriMacro$cell.line == "HCC1954"], "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("HCC1954")
plotReducedDim(sce.PeriMacro, "TSNE", colour_by="detected", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of detected genes")
plotReducedDim(sce.PeriMacro, "TSNE", colour_by="sum", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Nr. of total counts")
plotReducedDim(sce.PeriMacro, "TSNE", colour_by="Experiment", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Which experiment")
plotReducedDim(sce.PeriMacro[,sce.VECs$Experiment == "1"], "TSNE", colour_by="cyclone.phase", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Cycle stage in Experiment 1")
plotReducedDim(sce.PeriMacro[,sce.VECs$Experiment == "2"], "TSNE", colour_by="cyclone.phase", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("Cycle stage in Experiment 2")
plotReducedDim(sce.PeriMacro, "TSNE", colour_by="cell.line", text_colour = "red", text_by="subclusters",text_size = 4,add_legend = TRUE) + ggtitle("cell line")
dev.off()





#####
# bar plots VECs
#
#

cluster <- sce.PeriMacro$subclusters
condition <- sce.PeriMacro$cell.line
location <- str_remove(sce.PeriMacro$SampleName,"_.")
Experiment <- sce.PeriMacro$Experiment
data <- data.frame(cluster,condition,Experiment,location)
# Stacked

table(cluster,condition)
max.cells <- table(condition)
# HCC1954  MDA231 
#  503     278

max.cells.Exp <- table(Experiment)
# 1     2 
#601 180

data.norm <- data
data.norm$scaling.condition <-  data.norm$scaling.Experiment  <-  data.norm$scaling.location <- 1 

data.norm[data.norm$condition == "HCC1954",]$scaling.condition = 1/ max.cells["HCC1954"]
data.norm[data.norm$condition == "MDA231",]$scaling.condition = 1/ max.cells["MDA231"]

data.norm[data.norm$Experiment == 1,]$scaling.Experiment = 1/ max.cells.Exp[1]
data.norm[data.norm$Experiment == 2,]$scaling.Experiment = 1/ max.cells.Exp[2]




pdf("FIGURES/05_PeriMacro_Barplot_clusters.pdf", width = 5, height = 5)
ggplot(data.norm, aes(fill=location, y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, both cell lines")
ggplot(data.norm[data.norm$condition == "MDA231",], aes(fill=location , y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, Only MDA231") + scale_fill_manual(values=c("darkgreen",  "darkblue"))
ggplot(data.norm[data.norm$condition == "HCC1954",], aes(fill=location , y=scaling.condition, x=cluster)) +  
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, Only HCC1954") + scale_fill_manual(values=c("darkgreen",  "darkblue"))
ggplot(data.norm, aes(fill=condition, y=scaling.condition, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, cluster specificity wrt. cell lines")+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
ggplot(data.norm, aes(fill=Experiment, y=scaling.Experiment, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions, differences in Experiments")+ scale_fill_manual(values=c("lightgreen",  "#56B4E9"))
dev.off()


### save sce with fine annotation
saveRDS(sce,"Rdata/sce_fineAnnotation.rds")
