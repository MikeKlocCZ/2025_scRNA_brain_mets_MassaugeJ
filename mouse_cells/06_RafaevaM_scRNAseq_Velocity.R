
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
  library(bluster)
  library(velociraptor)
})

#====
sce.Microglia <- readRDS("Rdata/sceMicroglia_clusters.rds")

#clustering at the level of only Microglias
set.seed(101) 
nn.clusters.Microglia <- clusterCells(sce.Microglia, use.dimred="fastMNN",  BLUSPARAM=SNNGraphParam(k=25, type="rank", cluster.fun="walktrap"))
table(nn.clusters.Microglia)
#  1    2    3    4    5    6    7    8    9   10   11 
#1907  536  870 1836 1637 1337  117 1664  353 1475  327

#clsuter Microgila.2 is the one cycling
sce.Microglia$subclusters <- paste0("Microglia.",nn.clusters.Microglia)

#equip sce with the dimensioanl reduction
#sce <- runPCA(sce) #on log data!!
set.seed(101)
sce.Microglia <- runTSNE(sce.Microglia, dimred = "fastMNN")
set.seed(110)
sce.Microglia <- runUMAP(sce.Microglia, dimred =  "fastMNN", ncomponents=5, n_neighbors=25)

plotReducedDim(sce.Microglia, "TSNE", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)
plotReducedDim(sce.Microglia, "UMAP", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)



###########################
#trajectories
dec.Microglia <- modelGeneVar(sce.Microglia)
top.hvgs.Microglia <- getTopHVGs(dec.Microglia, n=2000)
velo.out.Microglia <- scvelo(sce.Microglia, subset.row=top.hvgs.Microglia, assay.X="spliced") #static model
velo.out.Microglia.dyn <- scvelo(sce.Microglia, subset.row=top.hvgs.Microglia, assay.X="spliced", mode = "dynamical") #dynamical model
saveRDS(velo.out.Microglia.dyn,"Rdata/dynVelocity_Microglia.rds")

sce.Microglia$velocity_pseudotime <- velo.out.Microglia$velocity_pseudotime
plotReducedDim(sce.Microglia,dim = "fastMNN", colour_by="velocity_pseudotime")

embedded.TSNE <- embedVelocity(reducedDim(sce.Microglia, "TSNE"), velo.out.Microglia)
grid.df <- gridVectors(sce.Microglia, embedded.TSNE, use.dimred = "TSNE")

embedded.UMAP <- embedVelocity(reducedDim(sce.Microglia, "UMAP"), velo.out.Microglia)
grid.df.UMAP <- gridVectors(sce.Microglia, embedded.UMAP, use.dimred = "UMAP")


pdf("plots/microglia_ssVelocityPlots.pdf", width = 10, height = 10)
plotReducedDim(sce.Microglia, "TSNE", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)
plotReducedDim(sce.Microglia, "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("both cell lines")

plotTSNE(sce.Microglia, colour_by="velocity_pseudotime") +
  geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, 
                                         xend=end.1, yend=end.2, colour=NULL), arrow=arrow(length=unit(0.05, "inches"))) 
plotReducedDim(sce.Microglia, "UMAP", colour_by="velocity_pseudotime") +
  geom_segment(data=grid.df.UMAP, mapping=aes(x=start.1, y=start.2, 
                                              xend=end.1, yend=end.2, colour=NULL), arrow=arrow(length=unit(0.01, "inches"))) 
dev.off()

sce.Microglia$dyn.velocity_pseudotime <- velo.out.Microglia.dyn$velocity_pseudotime
plotReducedDim(sce.Microglia,dim = "fastMNN", colour_by="dyn.velocity_pseudotime")

embedded.TSNE.dyn <- embedVelocity(reducedDim(sce.Microglia, "TSNE"), velo.out.Microglia.dyn)
grid.df.dyn <- gridVectors(sce.Microglia, embedded.TSNE.dyn, use.dimred = "TSNE")

embedded.UMAP.dyn <- embedVelocity(reducedDim(sce.Microglia, "UMAP"), velo.out.Microglia.dyn)
grid.df.UMAP.dyn <- gridVectors(sce.Microglia, embedded.UMAP.dyn, use.dimred = "UMAP")


pdf("plots/microglia_dynVelocityPlots.pdf", width = 10, height = 10)
plotReducedDim(sce.Microglia, "TSNE", colour_by="subclusters", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE)
plotReducedDim(sce.Microglia, "TSNE", colour_by="SampleName", text_colour = "black", text_by="subclusters",text_size = 4,add_legend = TRUE) +ggtitle("both cell lines")

plotTSNE(sce.Microglia, colour_by="dyn.velocity_pseudotime") +
  geom_segment(data=grid.df.dyn, mapping=aes(x=start.1, y=start.2, 
                                         xend=end.1, yend=end.2, colour=NULL), arrow=arrow(length=unit(0.05, "inches"))) 

plotReducedDim(sce.Microglia, "UMAP", colour_by="dyn.velocity_pseudotime") +
  geom_segment(data=grid.df.UMAP.dyn, mapping=aes(x=start.1, y=start.2, 
                                              xend=end.1, yend=end.2, colour=NULL), arrow=arrow(length=unit(0.01, "inches"))) 
dev.off()


#saveRDS(sce.Microglia,"Rdata/sceMicroglia_VelocityTests.rds")
sce.Microglia <- readRDS("Rdata/sceMicroglia_VelocityTests.rds")


#back to the top variable
top.hvgs.Microglia 

assay(velo.out.Microglia.dyn, "velocity")[1:20,1:4]
assay(velo.out.Microglia, "velocity")[1:20,1:4] #everything

head(dec.Microglia[top.hvgs.Microglia ,])

plotVelocity(velo.out.Microglia.dyn , c("Ccl4","Il1b","Btg2"))
plotVelocity(velo.out.Microglia , c("Ccl4","Il1b","Btg2"))

#assay(velo.out.Microglia.dyn, "velocity")[4,] %>% is.na %>% all

genes.all.nas <- apply(assay(velo.out.Microglia.dyn, "velocity"),1, function(ii) all(is.na(ii)))
genes.to.check <- top.hvgs.Microglia[top.hvgs.Microglia %in% names(genes.all.nas[which(genes.all.nas == FALSE)])]


pdf("plots/phase_plots.pdf", width = 10, height = 5)
for (gene in genes.to.check[1:150]) {
print(plotVelocity(velo.out.Microglia.dyn, gene))
}
dev.off()
