

##
#  QC on single-cell data (human part with cancer cells)
#  from Masague paper, GSE223309

#Use Bioconductor 3.21 (R 4.5)
if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(version = "3.21" ,force = TRUE)


pkg <- c("SingleCellExperiment",  "edgeR",
         "scater", "scran", "umap", "DropletUtils", "tidyverse", "cowplot",
         "BiocParallel", "BiocSingular")

pkg <- setdiff(pkg, rownames(installed.packages()))
if (length(pkg) > 0) {
  BiocManager::install(pkg)
}
## Load packages
suppressPackageStartupMessages({
  invisible(lapply(pkg, library, character.only = TRUE))
})


###################
#Make sure the directories exist, or create them
directories <- c("FIGURES","tables","Rdata")

for (directory in directories){
  if (!dir.exists(directory)){
    dir.create(directory)
    print(paste("dir", directory, "created"))
  } else {
    print(paste("dir", directory, "exists"))
  }
}



# Read-in data  MOUSE
sce <- readRDS("input_data/brainMets_cancer_sce.rds")
dim(sce )
#39197 18107

# STEP 1-> remove empty drops
table(colData(sce)$emptyDropsCR.isEmpty)
#FALSE  TRUE 
#17063  1044 

sce <- sce[,!colData(sce)$emptyDropsCR.isEmpty]
#  39197 17063

#Split By Experiment
colData(sce)$Experiment <-  str_remove(colData(sce)$SampleName, ".*_")



altExpNames(sce)
altExp(sce) #all 71 HTOs

#hashtag part, cut off those we cannot assign...investigate channel by channel



#used
hashtags <- c("HTO_A0255"="MDA231Br_cancer",
              "HTO_A0307"= "MDA231Br_stroma",
              "HTO_A0253" ="HCC1954Br_cancer",
              "HTO_A0306" ="HCC1954Br_stroma",
              "HTO_A0258" ="MDA231Br_cancer",
              "HTO_A0308" ="MDA231Br_stroma",
              "HTO_A0256" ="HCC1954Br_cancer")

altExp(sce) <- altExp(sce)[rownames(altExp(sce)) %in% names(hashtags),]



sce2 <- sce[,colData(sce)$Experiment == "2"]
sce1 <- sce[,colData(sce)$Experiment == "1"]


#
hto.counts.s1 <- as.matrix(counts(altExp(sce1)))
hto.counts.s1 <- hto.counts.s1[c("HTO_A0253","HTO_A0255"),]
rownames(hto.counts.s1) <- c("HCC1954","MDA231")

hto.counts.s2 <- as.matrix(counts(altExp(sce2)))
hto.counts.s2 <- hto.counts.s2[c("HTO_A0258","HTO_A0256"),]
rownames(hto.counts.s2) <- c("MDA231","HCC1954")


#logUMI vs logHTO
logUMIs <- log10(colSums(counts(sce1)) +1)
logHTOs <- log10(colSums(hto.counts.s1) +1)

df.plot <- data.frame(logUMIs = logUMIs,
                      logHTOs = logHTOs,
                      compartment = colData(sce1)$SampleName)

pdf("FIGURES/01_Exp1_HtosVsUmis.pdf")
ggplot(df.plot, aes(logUMIs,logHTOs,color = compartment) ) + geom_point() + ggtitle("Exp1, HTOs vs UMI exprs")
dev.off()

logUMIs.2 <- log10(colSums(counts(sce2)) +1)
logHTOs.2 <- log10(colSums(hto.counts.s2) +1)

df.plot.2 <- data.frame(logUMIs.2 = logUMIs.2,
                      logHTOs.2 = logHTOs.2,
                      compartment = colData(sce2)$SampleName)

pdf("FIGURES/01_Exp2_HtosVsUmis.pdf")
ggplot(df.plot.2, aes(logUMIs.2,logHTOs.2,color = compartment) ) + geom_point() + ggtitle("Exp2, HTOs vs UMI exprs")
dev.off()


#plot HTOs

df.plot.HTOs.1 <- data.frame(log10(t(hto.counts.s1)+1),
                             compartment = colData(sce1)$SampleName)


df.plot.HTOs.2 <- data.frame(log10(t(hto.counts.s2)+1),
                             compartment = colData(sce2)$SampleName)

pdf("FIGURES/01_HTOs_raw.pdf")
ggplot(df.plot.HTOs.1, aes(HCC1954,MDA231,color = compartment) ) + geom_point() + ggtitle("Exp1, HTOs")
ggplot(df.plot.HTOs.2, aes(HCC1954,MDA231,color = compartment) ) + geom_point() + ggtitle("Exp2, HTOs")
dev.off()



# BAD CELLS REMOVAL
#CHECK Mitochondrial content
#quality control: abundance of mitochondrial RNA + total counts
# Mitochondrial Reads


Mt <- rownames(sce)[which(str_detect(rowData(sce)$SYMBOL,"^MT-."))]
df <- perCellQCMetrics(sce,  subsets=list(Mito=Mt)) #exclude the tags

QC.lib <- isOutlier(df$sum, type="lower",nmads = 3) 
sum(QC.lib)
# 0 

QC.mit <- isOutlier(df$subsets_Mito_percent,  type="higher", nmads = 3) 
sum(QC.mit)
# 1340

attr(QC.mit,"thresholds")
#  lower  higher 
#-Inf 12.19359 

QC.expr <- isOutlier(df$detected, log=TRUE, type="lower",nmads = 3) 
sum(QC.expr)
# 911

attr(QC.expr,"thresholds")
#lower   higher 
#1864.445      Inf 

#cells that would be thrown away only for one reason 
sum(QC.mit != QC.expr)
#   1285

colData(sce) <- cbind(colData(sce), df)
sce$discard <- as.logical(QC.mit+QC.expr)

pdf("FIGURES/01_dicard_genes_CANCER.pdf", width = 4, height = 6)
par(mfrow=c(3,1))
gridExtra::grid.arrange(
  plotColData(sce, x="Experiment", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="Experiment", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="Experiment", y="subsets_Mito_percent", colour_by="discard") +  
    ggtitle("Mito percent"),
  ncol=1
) 
dev.off()

pdf("FIGURES/01_Mito_discard_CANCER.pdf", width = 8, height = 4)
par(mfrow=c(1,2))
plotColData(sce, x="sum", y="subsets_Mito_percent", 
            colour_by="discard", other_fields=c( "Experiment")) +
  facet_grid(~Experiment) +
  theme(panel.border = element_rect(color = "grey"))
dev.off()


#investigate what I am throwing away
sce.damaged <- sce[,sce$discard]
sce.damaged <- logNormCounts(sce.damaged)

lost <- calculateAverage(counts(sce.damaged))
kept <- calculateAverage(counts(sce[,!sce$discard]))


logged <- edgeR::cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

pdf("FIGURES/01_MA_discarded_cells_CANCER.pdf", width = 6, height = 4)
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[Mt], logFC[Mt], col="dodgerblue", pch=16)
text(abundance[Mt], logFC[Mt], labels = "Mito")
points(abundance["Malat1"], logFC["Malat1"], col="orange", pch=16)
text(abundance["Malat1"], logFC["Malat1"], labels = "Malat1")
dev.off()

sce <- sce[,!sce$discard]
dim(sce)
# 22542 15295

## Save filtered Data to my dir
#saveRDS(sce,"Rdata/sce_filtered.rds")



###
# AFter cleanup


sce2 <- sce[,colData(sce)$Experiment == "2"]
sce1 <- sce[,colData(sce)$Experiment == "1"]


#
hto.counts.s1 <- as.matrix(counts(altExp(sce1)))
hto.counts.s1 <- hto.counts.s1[c("HTO_A0253","HTO_A0255"),]
rownames(hto.counts.s1) <- c("HCC1954","MDA231")

hto.counts.s2 <- as.matrix(counts(altExp(sce2)))
hto.counts.s2 <- hto.counts.s2[c("HTO_A0258","HTO_A0256"),]
rownames(hto.counts.s2) <- c("MDA231","HCC1954")

#logUMI vs logHTO
logUMIs <- log10(colSums(counts(sce1)) +1)
logHTOs <- log10(colSums(hto.counts.s1) +1)

df.plot <- data.frame(logUMIs = logUMIs,
                      logHTOs = logHTOs,
                      compartment = colData(sce1)$SampleName)

pdf("FIGURES/01_Exp1_HtosVsUmis_afterQC.pdf")
ggplot(df.plot, aes(logUMIs,logHTOs,color = compartment) ) + geom_point() + ggtitle("Exp1, HTOs vs UMI exprs")
dev.off()

logUMIs.2 <- log10(colSums(counts(sce2)) +1)
logHTOs.2 <- log10(colSums(hto.counts.s2) +1)

df.plot.2 <- data.frame(logUMIs.2 = logUMIs.2,
                        logHTOs.2 = logHTOs.2,
                        compartment = colData(sce2)$SampleName)

pdf("FIGURES/01_Exp2_HtosVsUmis_afterQC.pdf")
ggplot(df.plot.2, aes(logUMIs.2,logHTOs.2,color = compartment) ) + geom_point()  + ggtitle("Exp2, HTOs vs UMI exprs")
dev.off()


#plot HTOs

df.plot.HTOs.1 <- data.frame(log10(t(hto.counts.s1)+1),
                             compartment = colData(sce1)$SampleName)


df.plot.HTOs.2 <- data.frame(log10(t(hto.counts.s2)+1),
                             compartment = colData(sce2)$SampleName)

pdf("FIGURES/01_HTOs_afterQC.pdf")
ggplot(df.plot.HTOs.1, aes(HCC1954,MDA231,color = compartment) ) + geom_point() + ggtitle("Exp1, HTOs")
ggplot(df.plot.HTOs.2, aes(HCC1954,MDA231,color = compartment) ) + geom_point() + ggtitle("Exp2, HTOs")
dev.off()



#remove those with zeros
hto.counts.s1 <- hto.counts.s1[,!(colSums(hto.counts.s1) == 0)]
hto.counts.s2 <- hto.counts.s2[,!(colSums(hto.counts.s2) == 0)]


dim(hto.counts.s1)
#     2 8901
dim(hto.counts.s2)
#    2 6394


summary(colSums(hto.counts.s1))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  43     538    2001    5582    7076  111630 
summary(colSums(hto.counts.s2))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#84    1754    5236    9856   12860  760711 

#demultiplexing + doublets base on HTOs
# rund directly hashedDrops...it will assign
#it establishes ambient profile of HTOs and then test the cells on that...should do it SampleSpecific
set.seed(101)
hash.stats.s1 <- hashedDrops(hto.counts.s1,constant.ambient=TRUE) #have only 2 conds
#hash.stats.s1b <- hashedDrops(hto.counts.s1,constant.ambient=FALSE)
set.seed(110)
hash.stats.s2 <- hashedDrops(hto.counts.s2,constant.ambient=TRUE)

df.plot.HTOs.1 <- df.plot.HTOs.1[colnames(hto.counts.s1),]
df.plot.HTOs.1$Confident <- hash.stats.s1$Confident

df.plot.HTOs.2 <- df.plot.HTOs.2[colnames(hto.counts.s2),]
df.plot.HTOs.2$Confident <- hash.stats.s2$Confident

pdf("FIGURES/01_HTOs_afterQC_confidence.pdf")
ggplot(df.plot.HTOs.1, aes(HCC1954,MDA231,color = Confident) ) + geom_point() + ggtitle("Exp1, HTOs")
ggplot(df.plot.HTOs.2, aes(HCC1954,MDA231,color = Confident) ) + geom_point() + ggtitle("Exp2, HTOs")
dev.off()

#check the assigned ambient profiles...DIFFERENT!
metadata(hash.stats.s1)$ambient
#  HCC1954   MDA231 
#210.9373 239.6180 

metadata(hash.stats.s2)$ambient
# MDA231   HCC1954 
#67.85629 240.45383 


c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,223, max = 255, alpha = 120, names = "lt.pink")

h1 <-  hist(log10(hto.counts.s1[1,]), breaks = 30)
h1line <- abline(v=log10(metadata(hash.stats.s1)$ambient[1]), col="red", lty=2)

h2 <-  hist(log10(hto.counts.s1[2,]), breaks = 30)
h2line <- abline(v=log10(metadata(hash.stats.s1)$ambient[2]), col="red", lty=2)

h3 <-  hist(log10(hto.counts.s2[2,]), breaks = 30)
h3line <- abline(v=log10(metadata(hash.stats.s2)$ambient[2]), col="red", lty=2)

h4 <-  hist(log10(hto.counts.s2[1,]), breaks = 30)
h4line <- abline(v=log10(metadata(hash.stats.s2)$ambient[1]), col="red", lty=2)


pdf("FIGURES/01_HTOs_profiles_ambient_estimation_CANCER.pdf")
par(mfrow=c(2,2))
gridExtra::grid.arrange(
  plot(h1, col = c1, xlab="Log[10] HCC1954 counts, Exp 1", main=""),
  abline(v=log10(metadata(hash.stats.s1)$ambient[1]), col="red", lty=2),
  plot(h2, col = c1, xlab="Log[10] MDA231 counts, Exp 1", main=""),
  abline(v=log10(metadata(hash.stats.s1)$ambient[2]), col="red", lty=2),
  plot(h3, col = c1, xlab="Log[10] HCC1954 counts, Exp 2", main=""),
  abline(v=log10(metadata(hash.stats.s2)$ambient[2]), col="red", lty=2),
  plot(h4, col = c1, xlab="Log[10] MDA231 counts, Exp 2", main=""),
  abline(v=log10(metadata(hash.stats.s2)$ambient[1]), col="red", lty=2),
  ncol=2
)
dev.off()

#how many cells we have?
table(hash.stats.s1$Confident) 
#FALSE  TRUE 
# 2711  6190 

table(hash.stats.s2$Confident) 
#FALSE  TRUE 
#537  5857 

h1.conf <-  hist(log10(hto.counts.s1[1,hash.stats.s1$Confident & hash.stats.s1$Best == 1]), breaks = 30)
h2.conf <-  hist(log10(hto.counts.s1[2,hash.stats.s1$Confident & hash.stats.s1$Best == 2]), breaks = 30)


pdf("FIGURES/01_HTOs_assignment_exp1_CANCER.pdf")
par(mfrow=c(1,2))
gridExtra::grid.arrange(
  plot(h1.conf, col = c1, xlab="Log[10] HCC1954 counts, Exp 1", main="Assigned to condition"),
  abline(v=log10(metadata(hash.stats.s1)$ambient[1]), col="red", lty=2),
  plot(h2.conf, col = c1, xlab="Log[10] MDA231 counts, Exp 1", main="Assigned to condition"),
  abline(v=log10(metadata(hash.stats.s1)$ambient[2]), col="red", lty=2),
  ncol=2
)
dev.off()



#mostly low count data, those which were not distinguishable with logFC1

hto.confident.s1 <- hto.counts.s1[,hash.stats.s1$Confident]
hto.confident.s2 <- hto.counts.s2[,hash.stats.s2$Confident]
summary(colSums(hto.confident.s1))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#60    1641    4062    7665   10479  111630
summary(colSums(hto.confident.s2))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#132    2251    5942   10510   13481  760711


r.s1 <- rank(-hash.stats.s1$Total)
color_Confident.s1 = c("black","green")[as.factor(hash.stats.s1$Confident )]

r.s2 <- rank(-hash.stats.s2$Total)
color_Confident.s2 = c("black","green")[as.factor(hash.stats.s2$Confident )]

pdf("FIGURES/01_HashTag_RankExp1_CANCER.pdf", width = 10, height = 5)
plot(r.s1, hash.stats.s1$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Confidently assigned, Exp1",col=color_Confident.s1 )
dev.off()
pdf("FIGURES/01_HashTag_RankExp2_CANCER.pdf", width = 10, height = 5)
plot(r.s2, hash.stats.s2$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Confidently assigned, Exp2",col=color_Confident.s2 )
dev.off()


#HereExp 1
#total hashcount over 100
#metadata(hash.stats.s1)$ambient %>% sum
hash.stats.s1$Confident[hash.stats.s1$Total < 450] %>% table
#FALSE  TRUE 
#1825   106  ...we do not trust these
hash.stats.s1$Confident[which(hash.stats.s1$Total < 430)] <- FALSE

#####################################################
#total cells Exp1
table(hash.stats.s1$Best[hash.stats.s1$Confident]) 
#  1    2 
#2722 3376

rownames(hto.counts.s1)
#1 = "HCC1954", 2 =  "MDA231"

# Hashtag Antibodies
Ht <- c("HCC1954", "MDA231")

#only confident:
hash.stats.s1 <- hash.stats.s1[hash.stats.s1$Confident,]


#IMPORTANT STEP: select only those cells confidently assigned to the wanted conditions
sce1 <- sce1[ , rownames(hash.stats.s1)]
sce1$cell.line <- Ht[hash.stats.s1$Best]


#Here Exp2
#total hashcount over 100
#metadata(hash.stats.s2)$ambient %>% sum
hash.stats.s2$Confident[hash.stats.s2$Total < 300] %>% table
#FALSE  TRUE 
#302    28   ...we do not trust these
hash.stats.s2$Confident[which(hash.stats.s2$Total < 300)] <- FALSE

#also from the HTO rank graph ... the first4 : doublets
#doublets.remove.s2 <- r.s2 %>% sort(decreasing = TRUE) %>% head(5) %>% names
#hash.stats.s2[doublets.remove.s2,]$Confident <- FALSE

#####################################################
#total cells Exp1
table(hash.stats.s2$Best[hash.stats.s2$Confident]) 
#    1    2 
#1897 3932 

rownames(hto.counts.s2)
#1 = "MDA231", 2 =  "HCC1954"

# Hashtag Antibodies
Ht.2 <- c("MDA231", "HCC1954")
#only confident:
hash.stats.s2 <- hash.stats.s2[hash.stats.s2$Confident,]


#IMPORTANT STEP: select only those cells confidently assigned to the wanted conditions
sce2 <- sce2[ , rownames(hash.stats.s2)]
sce2$cell.line <- Ht.2[hash.stats.s2$Best]





sce <- cbind(sce1,sce2) #Combine back


## Save Data to my dir
#saveRDS(sce,"Rdata/sce_demultiplexed_cleaned.rds")



