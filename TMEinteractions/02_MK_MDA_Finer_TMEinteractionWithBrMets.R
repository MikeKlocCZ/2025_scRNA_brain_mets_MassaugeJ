## Set top directory
#Sys.setenv(DISPLAY="localhost:15.0")


mydir <- "/scicore/home/bentires/GROUP/michal/rafaeva__maria/GSE223309_2023/human_data/starsolo_mapping_workflow/TMEinteractions"    # CHANGE!! working dir
setwd(mydir)

## Load packages
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidyverse)
  library(pheatmap)
  library(Seurat)
  library(scran)
  library(RColorBrewer)
  library(cowplot)
  library(ggpubr)
  library(scater)
  library(velociraptor)
})

#library(devtools)
#devtools::install_github("saeyslab/nichenetr",lib="/scicore/home/bentires/myberi81/R/x86_64-pc-linux-gnu-library/4.3")

#SETP 1: read the nichenet model
nichenet.dir <-  "/scicore/home/bentires/GROUP/michal/NicheNet_resources"   
lr_network <- readRDS(file.path(nichenet.dir,"human/lr_network.human.rds" ))
ligand_target_matrix <- readRDS(file.path(nichenet.dir,"human/ligand_target_matrix.human.rds" ))
weighted_networks <- readRDS(file.path(nichenet.dir,"human/weighted_networks.human.rds" ))


lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
##                  A2M        AANAT        ABCA1          ACE        ACE2
#A-GAMMA3'E 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
#1BG       0.0018503922 0.0011108718 0.0014225077 0.0028594037 0.001139013
#A1BG-AS1   0.0007400797 0.0004677614 0.0005193137 0.0007836698 0.000375007
#A1CF       0.0024799266 0.0013026348 0.0020420890 0.0047921048 0.003273375
#A2M        0.0084693452 0.0040689323 0.0064256379 0.0105191365 0.005719199

weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
#1 A-GAMMA3'E ACTG1P11   0.100 
#2 A-GAMMA3'E AXIN2      0.0869
#3 A-GAMMA3'E BUB1B-PAK6 0.0932
#4 A-GAMMA3'E CEACAM7    0.0793
#5 A-GAMMA3'E CHRNA1     0.0901
#6 A-GAMMA3'E DTX2P1     0.0976

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

library(nichenetr)
#: FROM stroma (mouse), TO cancer (human)
#convert genes to mouse (1-to-1 orthologs)
#lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = to) %>% drop_na()
#colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
#ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
#weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = to) %>% drop_na()


#
# get the relevant gene sets
#

## receiver ..CANCER
receiver <- "cancer"
expressed_genes_receiver <- read.csv("MDAbr_expressed_genes.csv",header = T)[,1] #just gene names, symbols
expressed_genes_receiver <- convert_alias_to_symbols(expressed_genes_receiver , "human", verbose = FALSE)
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] #3217


#diff exprs. genenes BrMet vs maternal
# MDA

#extract table from pdf
library(readxl)
geneset_oi  <- read_excel("DiffExp_MDA_BtVsParental.xlsx",col_names = FALSE) %>% unlist %>% unique #242
               
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)] 
geneset_oi <- convert_alias_to_symbols(geneset_oi, "human", verbose = FALSE)

## sender .. stroma
sender_celltypes <-  "TME"
dir.TME <- "/scicore/home/bentires/GROUP/michal/rafaeva__maria/GSE223309_2023/mouse_data/starsolo"
sce.TME <- readRDS(file.path(dir.TME,"Rdata/sce_fineAnnotation.rds")) #reduce to the cell types we are interested in (increased presence in in TME)

selection <- c("Vascular endothelial cells, capillary" , "Vascular endothelial cells, venous" , "Microglia", "Olfactory ensheathing cells",
               "Perivascular macrophages", "Oligodendrocytes precursor cells" , "Oligodendrocytes" , "Fibroblasts"  )
sce.TME <- sce.TME[,sce.TME$clusters %in% selection]

sce.TME[,sce.TME$clusters == "Vascular endothelial cells, capillary"]$FinerClusters %>% table
sce.TME[,sce.TME$clusters == "Vascular endothelial cells, venous"]$FinerClusters %>% table

sce.TME$clusters <- as.character(sce.TME$clusters)
sce.TME$clusters[sce.TME$clusters %in% c("Vascular endothelial cells, capillary" , "Vascular endothelial cells, venous")] <- "VECs"
sce.TME$clusters <- as.factor(sce.TME$clusters)



#cell.types <- c("VECs", "Microglia","Perivascular macrophages")
cell.types <- c("VECs", "Microglia","Perivascular macrophages")
for(cell.type in cell.types){
  sce.cell.type <- sce.TME[, sce.TME$clusters == cell.type]
  colLabels(sce.cell.type ) <- sce.cell.type$FinerClusters
  
  #filter unexpressed genes
  keep_feature <- rowSums(counts(sce.cell.type) > 5) > 10
  #
  sum(keep_feature) #10016
  sce.cell.type <- sce.cell.type[keep_feature,]
  rowData(sce.cell.type)$hugo.symbols <- rownames(sce.cell.type) %>% convert_mouse_to_human_symbols %>%
    convert_alias_to_symbols( "human", verbose = FALSE)
  
  
  expressed_genes_sender <- rowData(sce.cell.type)$hugo.symbols %>% na.omit %>% unique   # 2565
  #expressed_genes_sender <- convert_alias_to_symbols(expressed_genes_sender, "human", verbose = FALSE)
  
  sce.cell.type <- sce.cell.type[!is.na(rowData(sce.cell.type)$hugo.symbols),]
  
  condition_oi = "brain_mets"
  condition_reference = "parental" 
  
  ligands <- lr_network %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands,expressed_genes_sender)
  
  receptors <- lr_network %>% pull(to) %>% unique()
  expressed_receptors <- intersect(receptors,expressed_genes_receiver)
  
  lr_network_expressed <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
  #head(lr_network_expressed)
  
  #Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest
  potential_ligands <- lr_network_expressed %>% pull(from) %>% unique()
  
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities %>% arrange(-aupr_corrected) 
  
  best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
 # head(best_upstream_ligands)
  
  write.csv(ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected), paste0("top20_predictedLigandsMDA_",cell.type ,".csv"))
  
  # show histogram of ligand activity scores
  p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))), color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity ", y = "# ligands") +
    theme_classic()
  
  
  pdf(paste0("figures/MDA_ligands_histogram_",cell.type, "_on_mets.pdf"), width = 5, height = 5)
  p_hist_lig_activity + ggtitle(paste("MDA cancer cells affected by ",cell.type))
  dev.off()
  

  
  sce.sel <- sce.cell.type[rowData(sce.cell.type)$hugo.symbols %in% best_upstream_ligands,]
  sce.sel <- removeAltExps(sce.sel) 
  rownames(sce.sel) <- rowData(sce.sel)$hugo.symbols 
  
  #correct if more genes of the same names
  a <- table(rownames(sce.sel) )
  non.unique <- names(a[a>1])
  unique <- names(a[a==1])
  
  
  #keep the one with larger dispersion
  if (length(non.unique) > 0){
    sce.sel.unique <- sce.sel[unique,]
    sce.sel.nonunique <- sce.sel[ rowData(sce.sel)$hugo.symbols  %in% non.unique,]
    testStat <- apply(assays(sce.sel.nonunique)[["logcounts"]],1,IQR, na.rm = TRUE) # (inter-quantile distance) variance for each row
  #add rownames 
    names(testStat) <- paste(names(testStat),seq_along(testStat),sep="__")
    tSsp <- split.default(testStat, f = rownames(sce.sel.nonunique))
    kept <- sapply(tSsp, function(x) names(which.max(x)))
    head(kept) #if non unique ensemblids, keep those with higher variance
    ids.keep <- lapply(kept, str_remove,".*__") %>% unlist %>% as.numeric
    sce.sel.nonunique <-  sce.sel.nonunique[ids.keep,]#
  
    sce.sel <- rbind(sce.sel.unique,sce.sel.nonunique)}
  
  seuratObj <- as.Seurat(sce.sel,counts = "counts",
                         data = "logcounts")
  
  Idents(seuratObj)  <- colLabels(sce.sel)
  rotated.dotplot <-  DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
  pdf(paste0("figures/MDA_dotplot_",cell.type,".pdf"), width = 7, height = 4)
  print(rotated.dotplot )
  dev.off()
  
  
  #Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, 
                                                                   ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
  active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
  #nrow(active_ligand_target_links_df)
  #head(active_ligand_target_links_df)
  active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
  
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique()
  vis_ligand_target <- active_ligand_target_links[order_targets,order_ligands] %>% t()
  
  p_ligand_target_network <- vis_ligand_target %>% make_heatmap_ggplot(paste0("Prioritized ",cell.type, "-ligands"),"differentially expressed genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.1,0.2)) + theme(axis.text.x = element_text(face = "italic"))
  
  
  pdf(paste0("figures/MDA_target_regulatory_potential",cell.type, "_on_mets.pdf"), width = 8, height = 5)
  print(p_ligand_target_network)
  dev.off()
  
  
  ligand_aupr_matrix <- ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  vis_ligand_aupr <- ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr  <- vis_ligand_aupr %>% make_heatmap_ggplot(paste0("Prioritized ",cell.type,"-ligands"),"Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)")
  p_ligand_aupr
  
  # get the ligand-receptor network of the top-ranked ligands
  lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()
  
  # get the weights of the ligand-receptor interactions as used in the NicheNet model
  
  lr_network_top_df <- weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  # convert to a matrix
  lr_network_top_df <- lr_network_top_df %>% spread("from","weight",fill = 0)
  lr_network_top_matrix <- lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  # perform hierarchical clustering to order the ligands and receptors
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot(paste("Prioritized ",cell.type, "-ligands"),"Receptors expressed by malignant cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  
  
  pdf(paste0("figures/MDA_ligand_receptor_regulatory_potential",cell.type, "_on_mets.pdf"), width = 8, height = 5)
  print(p_ligand_receptor_network)
  dev.off()
  
  #link it to the individual gene expression
  
  #proportions:
  
  cluster <- colLabels(sce.cell.type)
  condition <- sce.cell.type$cell.line
  location <- str_remove(sce.cell.type$SampleName,"_.")
  data <- data.frame(cluster,condition,location)
  # Stacked
  
  
  max.cells <- table(condition)
  # HCC1954  MDA231 
  #x    y
  
  data.norm <- data
  data.norm$scaling.condition <-  1 
  
  data.norm[data.norm$condition == "HCC1954",]$scaling.condition = 1/ max.cells["HCC1954"]
  data.norm[data.norm$condition == "MDA231",]$scaling.condition = 1/ max.cells["MDA231"]
  
  barplot.location <- ggplot(data.norm[data.norm$condition == "MDA231",], aes(fill=location, y=scaling.condition, x=cluster)) + 
    geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Normalized proportions,MDA231")
  
  pdf(paste0("figures/MDA_Barplot_location_",cell.type,".pdf"), width = 8, height = 4.5)
   print(barplot.location)
  dev.off()
  
  #TSNES + velocity
  dec.cell.type <- modelGeneVar(sce.cell.type)
  top.hvgs.cell.type <- getTopHVGs(dec.cell.type, n=2000)
  velo.out.cell.type <- scvelo(sce.cell.type, subset.row=top.hvgs.cell.type, assay.X="spliced")
  
  sce.cell.type$velocity_pseudotime <- velo.out.cell.type$velocity_pseudotime
  embedded.TSNE <- embedVelocity(reducedDim(sce.cell.type, "TSNE"), velo.out.cell.type)
  grid.df <- gridVectors(sce.cell.type, embedded.TSNE, use.dimred = "TSNE")
  
  
  tsne1 <- plotReducedDim(sce.cell.type[,sce.cell.type$cell.line == "MDA231"], "TSNE", colour_by="FinerClusters", text_colour = "black", text_by="FinerClusters",text_size = 4,add_legend = FALSE) +ggtitle("MDA231")
  tsne.velo <- plotTSNE(sce.cell.type[,sce.cell.type$cell.line == "MDA231"], colour_by="velocity_pseudotime") +
    geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, 
                                           xend=end.1, yend=end.2, colour=NULL), arrow=arrow(length=unit(0.05, "inches"))) 
  
  #combined plot
  # Combine figures and legend separately
  figures_toprow <- cowplot::plot_grid(
    p_ligand_aupr + theme(legend.position = "right", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
   # p_ligand_aupr + theme(axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
     p_ligand_target_network + theme(legend.position = "bottom", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
    p_ligand_receptor_network + theme(legend.position = "bottom", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
    align = "hv",
    nrow = 1,
    rel_widths = c(60, 80, 80))
  
  #print(figures_toprow)
  
  figures_bottomprow <- cowplot::plot_grid(
    rotated.dotplot + theme(legend.position = "right", axis.ticks = element_blank(), axis.title.x = element_text(size = 12),
                            axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) +
      ylab("") + xlab("") + scale_y_discrete(position = "right"),
    barplot.location + theme( legend.position = "top", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
    align = "hv",
    nrow = 1,
    rel_widths = c(1, 0.8))
  #print(figures_bottomprow)
  
  figures_tsnes <- cowplot::plot_grid(
    tsne1 + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12),
                            axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) +
      ylab("") + xlab("") + scale_y_discrete(position = "right"),
    tsne.velo + theme( legend.position = "top", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
    align = "h",
    nrow = 1,
    rel_widths = c(1, 0.8))
  #print(figures_tsnes)
  
  combined_plot <- cowplot::plot_grid(figures_toprow, figures_bottomprow, figures_tsnes ,nrow = 3, align = "hv") 
  
  pdf(paste0("figures/MDA231_CombinedNicheNet_",cell.type,".pdf"), width = 18, height = 15)
  print(combined_plot)
  dev.off()
}

