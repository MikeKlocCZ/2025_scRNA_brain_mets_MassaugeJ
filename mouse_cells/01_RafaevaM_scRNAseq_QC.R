
##
#  QC on single-cell data 
#  from Masague paper


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
# STEP 1: upload the data from RDS format, after the STARsolo alignment


# Read-in data 
sce <- readRDS("results/join_sce/brainMets_mm10_ensdb102_sce.rds")

# STEP 2 -> assign condition  based on HTOs


altExpNames(sce)
altExp(sce)

#hashtag part, cut off those we cannot assign...investigate channel by channel

hto.counts <- as.matrix(counts(altExp(sce)))
dim(hto.counts)
#8 42844


#MOST DO NOT HAVE COUNTS
table(colSums(hto.counts) > 0, colData(sce)$SampleGroup)
#FALSE  TRUE 
#40204  2640 


#HTO_A0255 MDA231Br_cancer
#HTO_A0307 MDA231Br_stroma
#HTO_A0253 HCC1954Br_cancer
#HTO_A0306 HCC1954Br_stroma
#HTO_A0258 MDA231Br_cancer
#HTO_A0308 MDA231Br_stroma
#HTO_A0256 HCC1954Br_cancer
hashtags <- c("HTO_A0255"="MDA231Br_cancer",
             "HTO_A0307"= "MDA231Br_stroma",
             "HTO_A0253" ="HCC1954Br_cancer",
             "HTO_A0306" ="HCC1954Br_stroma",
             "HTO_A0258" ="MDA231Br_cancer",
             "HTO_A0308" ="MDA231Br_stroma",
             "HTO_A0256" ="HCC1954Br_cancer")
hto.counts <- hto.counts[names(hashtags),]
table(colSums(hto.counts) > 0, colData(sce)$SampleGroup)

#check scatter plots
df.plot <- as.data.frame(log2(t(hto.counts) +1))
df.plot <- df.plot |> mutate(Sample =  str_remove(rownames(df.plot),"-.*"))

pdf(paste0("plots/MDA231_HTOs_inMouseSamples.pdf"), width = 5,height = 5)
ggplot(df.plot, aes(x=HTO_A0255,y=HTO_A0307, color = Sample)) + geom_point() + xlab("anti_human HTO_A0255")  + ylab("anti_mouse HTO_A0307") +
  ggtitle("Exp1, MDA231 barcodes logexp")

ggplot(df.plot, aes(x=HTO_A0258,y=HTO_A0308, color = Sample)) + geom_point() + xlab("anti_human HTO_A0258")  + ylab("anti_mouse HTO_A0308") +
  ggtitle("Exp2, MDA231 barcodes logexp")
dev.off()

pdf(paste0("plots/HCC1954_HTOs_inMouseSamples.pdf"), width = 5,height = 5)
ggplot(df.plot, aes(x=HTO_A0253,y=HTO_A0306, color = Sample)) + geom_point() + xlab("anti_human HTO_A0253")  + ylab("anti_mouse HTO_A0306 (also used in Exp2)") +
  ggtitle("Exp1, HCC1954 barcodes logexp")

ggplot(df.plot, aes(x=HTO_A0256,y=HTO_A0306, color = Sample)) + geom_point() + xlab("anti_human HTO_A0256")  + ylab("anti_mouse HTO_A0308 (also used in Exp1)") +
  ggtitle("Exp2, HCC1954 barcodes logexp")
dev.off()

CDKN2A.raw.logcounts <- log2(assay(sce,"counts")["ENSMUSG00000044303",] + 1)
summary(CDKN2A.raw.logcounts)
df.plot$CDKN2A.deletion <- CDKN2A.raw.logcounts == 0
table(df.plot$CDKN2A.deletion)

pdf(paste0("plots/CDKN2AdeletionChecks_HTOs_inMouseSamples.pdf"), width = 5,height = 5)
ggplot(df.plot, aes(x=HTO_A0255,y=HTO_A0307, color = CDKN2A.deletion)) + geom_point(size = 1) + xlab("anti_human HTO_A0255")  + ylab("anti_mouse HTO_A0307") +
  ggtitle("Exp1, MDA231 barcodes logexp")
ggplot(df.plot, aes(x=HTO_A0258,y=HTO_A0308, color = CDKN2A.deletion)) + geom_point(size = 1) + xlab("anti_human HTO_A0258")  + ylab("anti_mouse HTO_A0308") +
  ggtitle("Exp2, MDA231 barcodes logexp")

ggplot(df.plot, aes(x=HTO_A0253,y=HTO_A0306, color = CDKN2A.deletion)) + geom_point(size = 1) + xlab("anti_human HTO_A0253")  + ylab("anti_mouse HTO_A0306 (also used in Exp2)") +
  ggtitle("Exp1, HCC1954 barcodes logexp")

ggplot(df.plot, aes(x=HTO_A0256,y=HTO_A0306, color = CDKN2A.deletion)) + geom_point(size = 1) + xlab("anti_human HTO_A0256")  + ylab("anti_mouse HTO_A0308 (also used in Exp1)") +
  ggtitle("Exp2, HCC1954 barcodes logexp")
dev.off()


#my demultiplexing
raw <- "results/STARsolo/Stroma_1/brainMets_mm10_ensdb102_sce.rds"
sce.name = f"results/{config['remove_empty_droplets']}/{{sample}}/{{output_tag}}_sce.rds"



#human
sce.human <- readRDS("/scicore/home/bentires/GROUP/michal/rafaeva__maria/GSE223309_2023/human_data/starsolo_mapping_workflow/results/join_sce/brainMets_hg38_ensdb104_sce.rds")
hto.counts.human <- as.matrix(counts(altExp(sce.human)))

table(colSums(hto.counts.human) > 0, colData(sce.human)$SampleGroup)

hto.counts.human <- hto.counts.human[names(hashtags),]
table(colSums(hto.counts.human) > 0, colData(sce.human)$SampleGroup)
df.hum.plot <- as.data.frame(log2(t(hto.counts.human+1)))
df.hum.plot <- df.hum.plot |> mutate(Sample =  str_remove(rownames(df.hum.plot),"-.*"))

pdf(paste0("plots/MDA231_HTOs_inHumanSamples.pdf"), width = 5,height = 5)
ggplot(df.hum.plot, aes(x=HTO_A0255,y=HTO_A0307, color = Sample)) + geom_point() + xlab("anti_human HTO_A0255")  + ylab("anti_mouse HTO_A0307") +
  ggtitle("Exp1, MDA231 barcodes logexp")

ggplot(df.hum.plot, aes(x=HTO_A0258,y=HTO_A0308, color = Sample)) + geom_point() + xlab("anti_human HTO_A0258")  + ylab("anti_mouse HTO_A0308") +
  ggtitle("Exp2, MDA231 barcodes logexp")
dev.off()

pdf(paste0("plots/HCC1954_HTOs_inHumanSamples.pdf"), width = 5,height = 5)
ggplot(df.hum.plot, aes(x=HTO_A0253,y=HTO_A0306, color = Sample)) + geom_point() + xlab("anti_human HTO_A0253")  + ylab("anti_mouse HTO_A0306 (also used in Exp2)") +
  ggtitle("Exp1, HCC1954 barcodes logexp")

ggplot(df.hum.plot, aes(x=HTO_A0256,y=HTO_A0306, color = Sample)) + geom_point() + xlab("anti_human HTO_A0256")  + ylab("anti_mouse HTO_A0308 (also used in Exp1)") +
  ggtitle("Exp2, HCC1954 barcodes logexp")
dev.off()


CDKN2A.raw.logcounts <- log2(assay(sce.human,"counts")["ENSG00000147889",] + 1)
summary(CDKN2A.raw.logcounts)
df.hum.plot$CDKN2A.deletion <- CDKN2A.raw.logcounts == 0

pdf(paste0("plots/CDKN2AdeletionChecks_HTOs_inHumanSamples.pdf"), width = 5,height = 5)
ggplot(df.hum.plot, aes(x=HTO_A0255,y=HTO_A0307, color = CDKN2A.deletion)) + geom_point() + xlab("anti_human HTO_A0255")  + ylab("anti_mouse HTO_A0307") +
  ggtitle("Exp1, MDA231 barcodes logexp")
ggplot(df.hum.plot, aes(x=HTO_A0258,y=HTO_A0308, color = CDKN2A.deletion)) + geom_point() + xlab("anti_human HTO_A0258")  + ylab("anti_mouse HTO_A0308") +
  ggtitle("Exp2, MDA231 barcodes logexp")

ggplot(df.hum.plot, aes(x=HTO_A0253,y=HTO_A0306, color = CDKN2A.deletion)) + geom_point() + xlab("anti_human HTO_A0253")  + ylab("anti_mouse HTO_A0306 (also used in Exp2)") +
  ggtitle("Exp1, HCC1954 barcodes logexp")

ggplot(df.hum.plot, aes(x=HTO_A0256,y=HTO_A0306, color = CDKN2A.deletion)) + geom_point() + xlab("anti_human HTO_A0256")  + ylab("anti_mouse HTO_A0308 (also used in Exp1)") +
  ggtitle("Exp2, HCC1954 barcodes logexp")
dev.off()


#let's trust the  mouse cells are really mouse (looks like it's ok)

# Read-in data  MOUSE
sce <- readRDS("results/join_sce/brainMets_mm10_ensdb102_sce.rds")

# STEP 2 -> assign condition  based on HTOs


altExpNames(sce)
altExp(sce) #all 71 HTOs

#hashtag part, cut off those we cannot assign...investigate channel by channeö



#used
hashtags <- c("HTO_A0255"="MDA231Br_cancer",
              "HTO_A0307"= "MDA231Br_stroma",
              "HTO_A0253" ="HCC1954Br_cancer",
              "HTO_A0306" ="HCC1954Br_stroma",
              "HTO_A0258" ="MDA231Br_cancer",
              "HTO_A0308" ="MDA231Br_stroma",
              "HTO_A0256" ="HCC1954Br_cancer")

altExp(sce) <- altExp(sce)[rownames(altExp(sce)) %in% names(hashtags),]


#for demultiplexing, use only the mouse HTOs
hto.counts <- counts(altExp(sce))[c("HTO_A0307","HTO_A0306","HTO_A0308"),]

#demultiplexing + doublets base on HTOs
# rund directly hashedDrops...it will assign
#it establishes ambient profile of HTOs and then test the cells on that...should do it SampleSpecific
set.seed(132)
hash.stats <- hashedDrops(hto.counts)


#check the assigned ambient profiles...DIFFERENT!
metadata(hash.stats)$ambient
# HTO_A0307  HTO_A0306  HTO_A0308 
#5.4135039 35.4890167  0.1554298 


c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")


h1 <-  hist(log10(hto.counts[1,]), breaks = 30)
h1line <- abline(v=log10(metadata(hash.stats)$ambient[1]), col="red", lty=2)

h2 <-  hist(log10(hto.counts[2,]), breaks = 30)
h2line <- abline(v=log10(metadata(hash.stats)$ambient[2]), col="red", lty=2)

h3 <-  hist(log10(hto.counts[3,]), breaks = 30)
h3line <- abline(v=log10(metadata(hash.stats)$ambient[3]), col="red", lty=2)


pdf("plots/HTOs_profiles_ambient_estimation.pdf",width = 8,height = 3)
par(mfrow=c(1,3))
gridExtra::grid.arrange(
  plot(h1, col = c1, xlab="Log[10] HT1 counts", main=""),
  abline(v=log10(metadata(hash.stats)$ambient[1]), col="red", lty=2),
  plot(h2, col = c1, xlab="Log[10] HT2 counts", main=""),
  abline(v=log10(metadata(hash.stats)$ambient[2]), col="red", lty=2),
  plot(h3, col = c1, xlab="Log[10] HT3 counts", main=""),
  abline(v=log10(metadata(hash.stats)$ambient[3]), col="red", lty=2),
  ncol=3
)
dev.off()

#how many cells we have?
sum(hash.stats$Confident) 
#20713

#mostly low count data, those which were not distinguishable with logFC1
hash.stats$LowCount <- LowCount <- hash.stats$Doublet == 0 & hash.stats$Confident == 0


#hard-coded cut-off of low HTO linked to the ambient...
#ambient is around 20 for each channel, say minimal count in total is 40
hard.cutoff.s1 <- hash.stats.s1[hash.stats.s1$Confident == TRUE & hash.stats.s1$Total < 40,]
hard.cutoff.s2 <- hash.stats.s2[hash.stats.s2$Confident == TRUE & hash.stats.s2$Total < 40,]
#add them to lowcounts as well
hash.stats.s1[rownames(hard.cutoff.s1),]$Confident <- FALSE
hash.stats.s2[rownames(hard.cutoff.s2),]$Confident <- FALSE
hash.stats.s1[rownames(hard.cutoff.s1),]$LowCount <- TRUE
hash.stats.s2[rownames(hard.cutoff.s2),]$LowCount <-TRUE 

#how many cells we have?
sum(hash.stats.s1$Confident) + sum(hash.stats.s2$Confident)
#13084

##########
pdf("plots/Demultiplexing.pdf", width = 8, height = 6)
par(mfrow=c(2,2))
r.s1 <- rank(-hash.stats.s1$Total)
r.s2 <- rank(-hash.stats.s2$Total)
color_doublets.s1 = c("black","red")[as.factor(hash.stats.s1$Doublet )]
color_doublets.s2 = c("black","red")[as.factor(hash.stats.s2$Doublet )]
color_low.s1 = c("black","red")[as.factor(hash.stats.s1$LowCount )]
color_low.s2 = c("black","red")[as.factor(hash.stats.s2$LowCount )]
plot(r.s1, hash.stats.s1$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Doublets from FC2",col=color_doublets.s1)
plot(r.s2, hash.stats.s2$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Doublets from FC2",col=color_doublets.s2)
plot(r.s1, hash.stats.s1$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Not confident",col=color_low.s1)
plot(r.s2, hash.stats.s2$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Not confident",col=color_low.s2)
dev.off()


#hist(log10(hash.stats.s1$Total), xlab="Log[10] HTO count", main="Sample 1")
#hist(log10(hash.stats.s2$Total), xlab="Log[10] HTO count", main="Sample 2")
#hist(hash.stats.s1$LogFC, xlab="Log fold-change from best to second HTO", main="")

#############
###some checking of the values
hto.low.s1 <- hto.counts.s1[,LowCount.s1]
hto.low.s2 <- hto.counts.s2[,LowCount.s2]

summary(colSums(hto.low.s1))
summary(colSums(hto.low.s2))

sortvec.s1 <- order(colSums(hto.low.s1))
sortvec.s2 <- order(colSums(hto.low.s2))

hto.low.s1 <- hto.low.s1[,sortvec.s1]
hto.low.s2 <- hto.low.s2[,sortvec.s2]

hto.doublets.s1 <- hto.counts.s1[,hash.stats.s1$Doublet] 
hto.doublets.s2 <- hto.counts.s2[,hash.stats.s2$Doublet] 

hto.confident.s1 <- hto.counts.s1[,hash.stats.s1$Confident]
hto.confident.s2 <- hto.counts.s2[,hash.stats.s2$Confident]
summary(colSums(hto.confident.s1))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#40.0   389.8   509.0   662.0   744.0 20439.0 
summary(colSums(hto.confident.s2))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#44.0   428.0   556.0   727.7   827.0 14497.0 
sort.conf.s1 <- order(colSums(hto.confident.s1))
sort.conf.s2 <- order(colSums(hto.confident.s2))

hto.confident.s1 <- hto.confident.s1[,sort.conf.s1]
hto.confident.s2 <- hto.confident.s2[,sort.conf.s2]
#####################################################
#total cells
table(hash.stats.s1$Best[hash.stats.s1$Confident]) + table(hash.stats.s2$Best[hash.stats.s2$Confident])
# 1    2    3    4  ....Air = 1, Smoke = 2
#3133 3423 3288 3240

table(Ht[hash.stats.s1$Best[hash.stats.s1$Confident]]) + table(Ht[hash.stats.s2$Best[hash.stats.s2$Confident]])
#Air   AirTumor      Smoke SmokeTumor 
#3133       3288       3423       3240 

#only confident:
hash.stats.s1 <- hash.stats.s1[hash.stats.s1$Confident,]
hash.stats.s2 <- hash.stats.s2[hash.stats.s2$Confident,]

hash.stats.s1 <- hash.stats.s1[hash.stats.s1$Best == 1  | hash.stats.s1$Best == 2 ,]
hash.stats.s2 <- hash.stats.s2[hash.stats.s2$Best == 1  | hash.stats.s2$Best == 2 ,]
table(Ht[hash.stats.s1$Best]) + table(Ht[hash.stats.s2$Best])
#Air Smoke 
#3133  3423

#Total nr of cells
sum(table(Ht[hash.stats.s1$Best]) + table(Ht[hash.stats.s2$Best]))
#6556
sum(table(Ht[hash.stats.s1$Best])) 

#IMPORTANT STEP: select only those cells confidently assigned to the wanted conditions
sce1 <- sce1[ , rownames(hash.stats.s1)]
table(colData(sce1)$dataset)
#Sample1 Sample2 
#3285       0 

sce2 <- sce2[ , rownames(hash.stats.s2)]
table(colData(sce2)$dataset)
#Sample1 Sample2 
#0    3271 
sce1$cell.cond <- Ht[hash.stats.s1$Best]
sce2$cell.cond <- Ht[hash.stats.s2$Best]


sce <- cbind(sce1,sce2)
#omit HTOs part
altExp(sce) <- altExp(sce)[Ab,]
#dim: 31053 6558 

## Save Data to my dir
#saveRDS(sce,"Rdata/sce_demultiplexed.rds")

# Read-in data 
#sce <- readRDS("Rdata/sce_demultiplexed.rds")

#CHECK Mitochondrial content
#quality control: abundance of mitochondrial RNA + total counts
# Mitochondrial Reads
Mt <- rownames(sce)[which(as.character(rowData(sce)$chromosome_name) == "chrM")]
df <- perCellQCMetrics(sce,  subsets=list(Mito=Mt)) #exclude the tags

colnames(df) <- gsub(" ", "_", colnames(df)) #aesthetic: don't want empty slots in strings

QC.lib <- isOutlier(df$sum, type="lower",nmads = 3) 
sum(QC.lib)
# 0 

QC.mit <- isOutlier(df$subsets_Mito_percent,  type="higher", nmads = 3) 
sum(QC.mit)
#375

attr(QC.mit,"thresholds")

QC.expr <- isOutlier(df$detected, log=TRUE, type="lower",nmads = 3) 
sum(QC.expr)
# 35 cells with low expression level
attr(QC.expr,"thresholds")
#259 genes is the lower bound

#cells that would be thrown away only for one reason 
sum(QC.mit != QC.expr)
#340
colData(sce) <- cbind(colData(sce), df)
sce$discard <- as.logical(QC.mit+QC.expr)

pdf("plots/dicard_genes.pdf", width = 4, height = 6)
par(mfrow=c(3,1))
gridExtra::grid.arrange(
  plotColData(sce, x="cell.cond", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="cell.cond", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="cell.cond", y="subsets_Mito_percent", colour_by="discard") +
    ggtitle("Mito percent"),
  ncol=1
)
dev.off()

pdf("plots/Mito_discard.pdf", width = 8, height = 4)
par(mfrow=c(1,2))
plotColData(sce, x="sum", y="subsets_Mito_percent", 
            colour_by="discard", other_fields=c( "cell.cond")) +
  facet_grid(~cell.cond) +
  theme(panel.border = element_rect(color = "grey"))
dev.off()


#investigate what I am throwing away
sce.damaged <- sce[,sce$discard]
sce.damaged <- logNormCounts(sce.damaged)

lost <- calculateAverage(counts(sce.damaged))
kept <- calculateAverage(counts(sce[,!sce$discard]))

library(edgeR)
logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

pdf("plots/MA_discarded_cells.pdf", width = 6, height = 4)
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[Mt], logFC[Mt], col="dodgerblue", pch=16)
text(abundance[Mt], logFC[Mt], labels = "Mito")
points(abundance["Malat1"], logFC["Malat1"], col="orange", pch=16)
text(abundance["Malat1"], logFC["Malat1"], labels = "Malat1")
dev.off()

MeanAbundance.damaged <- rowMeans(logcounts(sce.damaged))
topGenes.damaged.id <- order(MeanAbundance.damaged,decreasing = TRUE)[1:30]
topGenes.damaged <- rownames(sce.damaged)[topGenes.damaged.id]
MeanAbundance.damaged[topGenes.damaged.id]
MeanAbundance.damaged[topGenes.damaged.id]
#Malat1   mt-Co3   mt-Co1  mt-Cytb  mt-Atp6   mt-Co2   mt-Nd4  Gm42418 
#6.998473 5.582174 5.520293 5.277102 5.192514 5.028968 4.628490 4.412114 
#mt-Nd1   mt-Nd2  mt-Nd4l   mt-Nd5   S100a9     Xist     Actb   S100a8 
#4.200398 3.831379 2.756435 2.672404 2.226927 1.975786 1.928378 1.878631 

malat.data <-sce.damaged[rownames(sce.damaged)=="Malat1",]
hist(logcounts(malat.data)[1,], xlab="Log[10] malat1 count", main="", breaks = 30)


#Malat1 is frequently detected in poly-A captured RNA-seq data, 
#independent of protocol. The level of detection appears to be 
#cell type-specific and also shows some inverse correlation with 
#cell health. We have observed that dead/dying cells have higher 
#expression of Malat1.
sce <- sce[,!sce$discard]
## Save filtered Data to my dir
#saveRDS(sce,"Rdata/sce_filtered.rds")

#sce <- readRDS("Rdata/sce_filtered.rds")
table(sce$cell.cond)
#Air Smoke 
#2923  3258 