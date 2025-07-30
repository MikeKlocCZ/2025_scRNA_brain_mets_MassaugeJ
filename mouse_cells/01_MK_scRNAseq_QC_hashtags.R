
##
#  QC on single-cell data 
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

########
# Read-in data 
# These are mouse cells (TME and stroma) from the brain
# 10X data aligned with STARsolo, mm10 , ensdb102
sce <- readRDS("input_data/brainMets_sce.rds")

########
#STEP 1: Check the sample annotation
unique(colData(sce)$SampleName)
#"Stroma_2" "Stroma_1" "TME_2"    "TME_1" 

#TME: cells in proximity to the metastatic lesions
#Stroma: Distant cells

#The experiment was performed in 2 batches (underscore)

#split by the batch
colData(sce)$Experiment <-  str_remove(colData(sce)$SampleName, ".*_")

########
#STEP 3: Hashtags in altExp
# Hashtags were used to multiplex both breast cancer cell lines 
# (triple-neg brain-tropic  MDAMB231 and HER2+ brain-tropic HCC1954)

altExpNames(sce)
altExp(sce) #all 71 TotalSeq HTOs

#Hashing scheme (see the Methods)
# Experiment (Batch) - Model - Hashtag 
# 1 MDA231-BrM        TotalSeq-A0255  anti-human Hashtag 5 
#                     TotalSeq-A0307  anti-mouse Hashtag 7 
#  HCC1954-BrM        TotalSeq-A0253  anti-human Hashtag 3 
#                     TotalSeq-A0306  anti-mouse Hashtag 6
# 2 MDA231-BrM.       TotalSeq-A0258  anti-human Hashtag 8 
#                     TotalSeq-A0308  anti-mouse Hashtag 8 
#   HCC1954-BrM       TotalSeq-A0256  anti-human Hashtag 6 
#                     TotalSeq-A0306  anti-mouse Hashtag 6

# strategy (from the Methods):
# single-cell suspension of each brain metastasis sample was incubated with an anti-human Hashtag antibody (binding
#to human cancer cells) and an anti-mouse Hashtag antibody (binding to mouse cells),
# The Calcein Violet+ Annexin V+ live cells of
# each sample were FACS-sorted into three groups, i.e., the labeling GFP+ mCherry+ cancer cells, labeled GFP+ mCherry+ TME cells,
# and rest unlabeled GFP+ mCherry+ cells. The single-cell suspension of uninoculated mouse brain, which only contained GFP 
# mCherry+ cells, was used to set the gates for GFP+ or mCherry+ cells.

#DICTIONARY (from the paper)
hashtags <- c("HTO_A0255"="MDA231Br_cancer",
              "HTO_A0307"= "MDA231Br_stroma",
              "HTO_A0253" ="HCC1954Br_cancer",
              "HTO_A0306" ="HCC1954Br_stroma",
              "HTO_A0258" ="MDA231Br_cancer",
              "HTO_A0308" ="MDA231Br_stroma",
              "HTO_A0256" ="HCC1954Br_cancer")

altExp(sce) <- altExp(sce)[rownames(altExp(sce)) %in% names(hashtags),]

hto.counts <- as.matrix(counts(altExp(sce)))
dim(hto.counts)
#  7 42844

#check specificity of the HTOs in the scatter plots
df.plot <- as.data.frame(log2(t(hto.counts) +1))
df.plot <- df.plot |> mutate(Sample =  str_remove(rownames(df.plot),"-.*"))

pdf(paste0("FIGURES/01_MDA231_HTOs_inMouseSamples.pdf"), width = 5,height = 5)
ggplot(df.plot, aes(x=HTO_A0255,y=HTO_A0307, color = Sample)) + geom_point() + xlab("anti_human HTO_A0255")  + ylab("anti_mouse HTO_A0307") +
  ggtitle("Exp1, MDA231 barcodes logexp")

ggplot(df.plot, aes(x=HTO_A0258,y=HTO_A0308, color = Sample)) + geom_point() + xlab("anti_human HTO_A0258")  + ylab("anti_mouse HTO_A0308") +
  ggtitle("Exp2, MDA231 barcodes logexp")
dev.off()

pdf(paste0("FIGURES/01_HCC1954_HTOs_inMouseSamples.pdf"), width = 5,height = 5)
ggplot(df.plot, aes(x=HTO_A0253,y=HTO_A0306, color = Sample)) + geom_point() + xlab("anti_human HTO_A0253")  + ylab("anti_mouse HTO_A0306 (also used in Exp2)") +
  ggtitle("Exp1, HCC1954 barcodes logexp")

ggplot(df.plot, aes(x=HTO_A0256,y=HTO_A0306, color = Sample)) + geom_point() + xlab("anti_human HTO_A0256")  + ylab("anti_mouse HTO_A0308 (also used in Exp1)") +
  ggtitle("Exp2, HCC1954 barcodes logexp")
dev.off()

# APPARENLTY THE HTOs were not human or mouse specific
# Here we are dealing with mouse cells only -> anti-human HTOs should not be present
#luckily, the human cancer cells vs the murine cells were also fluorescently tagged
# red only -> human cells, red + green -> TME, unlabelled -> stroma  


#for demultiplexing, use only the mouse HTOs
# split by the experiments

#split by batches
sce2 <- sce[,colData(sce)$Experiment == "2"]
sce1 <- sce[,colData(sce)$Experiment == "1"]


#
hto.counts.s1 <- as.matrix(counts(altExp(sce1)))
hto.counts.s1 <- hto.counts.s1[c("HTO_A0306","HTO_A0307"),]
rownames(hto.counts.s1) <- c("HCC1954","MDA231")

hto.counts.s2 <- as.matrix(counts(altExp(sce2)))
hto.counts.s2 <- hto.counts.s2[c("HTO_A0308","HTO_A0306"),]
rownames(hto.counts.s2) <- c("MDA231","HCC1954")


# QC plot
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






#demultiplexing + doublets base on HTOs
# rund directly hashedDrops...it will assign
#it establishes ambient profile of HTOs and then test the cells on that...should do it SampleSpecific
set.seed(132)
hash.stats.s1 <- hashedDrops(hto.counts.s1 , constant.ambient = TRUE) #only 2 HTOs
hash.stats.s2 <- hashedDrops(hto.counts.s2 , constant.ambient = TRUE)

#name the Barcodes correctly
hash.stats.s1$cell.line <- ifelse(hash.stats.s1$Best == 1, "HCC1954", "MDA231" )
hash.stats.s2$cell.line <- ifelse(hash.stats.s2$Best == 2, "HCC1954", "MDA231" ) #here they are swapped


#check the assigned ambient profiles...DIFFERENT!
metadata(hash.stats.s1)$ambient
#  HCC1954   MDA231 
#38.58537 22.44611 

metadata(hash.stats.s2)$ambient
# MDA231  HCC1954 
# 15.51681 30.75386 


#########
# Define color
c1 <- rgb(173, 216, 230, max = 255, alpha = 80, names = "lt.blue")

# Open PDF device
pdf("FIGURES/01_HTOs_profiles_ambient_estimation.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))  # 2x2 layout

hist(log10(hto.counts.s1[1, ]), breaks = 30, col = c1,
     main = "Exp. 1, HCC1954", xlab = "Log10 HTO1 counts")
abline(v = log10(metadata(hash.stats.s1)$ambient[1]), col = "red", lty = 2)

hist(log10(hto.counts.s1[2, ]), breaks = 30, col = c1,
     main = "Exp. 1, MDA231", xlab = "Log10 HTO2 counts")
abline(v = log10(metadata(hash.stats.s1)$ambient[2]), col = "red", lty = 2)

hist(log10(hto.counts.s2[2, ]), breaks = 30, col = c1,
     main = "Exp. 2, HCC1954", xlab = "Log10 HTO2 counts")
abline(v = log10(metadata(hash.stats.s2)$ambient[2]), col = "red", lty = 2)

hist(log10(hto.counts.s2[1, ]), breaks = 30, col = c1,
     main = "Exp. 2, MDA231", xlab = "Log10 HTO1 counts")
abline(v = log10(metadata(hash.stats.s2)$ambient[1]), col = "red", lty = 2)

dev.off()


#how many cells we have?
sum(hash.stats.s1$Confident) + sum(hash.stats.s2$Confident)
#21698


##########
pdf("FIGURES/01_Demultiplexing.pdf", width = 8, height = 6)
par(mfrow=c(1,2))
r.s1 <- rank(-hash.stats.s1$Total)
r.s2 <- rank(-hash.stats.s2$Total)
color_low.s1 = c("black","red")[as.factor(hash.stats.s1$Confident )]
color_low.s2 = c("black","red")[as.factor(hash.stats.s2$Confident )]
plot(r.s1, hash.stats.s1$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Confident (exp1)",col=color_low.s1)
plot(r.s2, hash.stats.s2$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Confident (exp2)",col=color_low.s2)
dev.off()

#In exp. 1 some of he most expressed HTOs are not cofident --> quite likely doublets
# association is based on logFC between the 1st and 2nd most abundant barcode, in these it's very low

sel <- hash.stats.s1[hash.stats.s1$Confident == FALSE,]
sel <- sel[order(-sel$Total),]
head(sel)

# DataFrame with 6 rows and 7 columns
# Total      Best    Second      LogFC    LogFC2   Doublet Confident
# <numeric> <integer> <integer>  <numeric> <numeric> <logical> <logical>
#   TME_1-GGTTAACCAGCTCTGG        53148         2        NA 0.39770316        NA        NA     FALSE
# Stroma_1-GTTGTGAGTGAGATCG     45034         1        NA 1.23788534        NA        NA     FALSE
# Stroma_1-TTGGGCGGTACGGGAT     35017         2        NA 0.00436291        NA        NA     FALSE
# Stroma_1-AAAGAACAGAGGGCGA     30360         2        NA 0.28735066        NA        NA     FALSE
# Stroma_1-GTTGTAGTCTGGGCAC     25645         2        NA 1.88308726        NA        NA     FALSE
# Stroma_1-TGCGGGTCAAGCGGAT     21880         2        NA 0.87160746        NA        NA     FALSE

#only confident:
hash.stats.s1 <- hash.stats.s1[hash.stats.s1$Confident,]
hash.stats.s2 <- hash.stats.s2[hash.stats.s2$Confident,]


#IMPORTANT STEP: select only those cells confidently assigned to the wanted conditions
sce1 <- sce1[ , rownames(hash.stats.s1)]
sce2 <- sce2[ , rownames(hash.stats.s2)]

sce1$cell.line <- hash.stats.s1$cell.line
sce2$cell.line  <- hash.stats.s2$cell.line

######
# AT THIS STAGE WE HAVE DEMULTIPLEXED THE CELL LINES

sce <- cbind(sce1,sce2)


########
#STEP 4:  BAD CELLS REMOVAL
#
#CHECK Mitochondrial content
#quality control: abundance of mitochondrial RNA + total counts
# Mitochondrial Reads

Mt <- rownames(sce)[which(str_detect(rowData(sce)$SYMBOL,"mt-."))]
df <- perCellQCMetrics(sce,  subsets=list(Mito=Mt)) #exclude the tags

QC.lib <- isOutlier(df$sum, type="lower",nmads = 3) 
sum(QC.lib)
# 0 

QC.mit <- isOutlier(df$subsets_Mito_percent,  type="higher", nmads = 3) 
sum(QC.mit)
# 1868

attr(QC.mit,"thresholds")
#  lower  higher 
#-Inf 8.718131 

QC.expr <- isOutlier(df$detected, log=TRUE, type="lower",nmads = 3) 
sum(QC.expr)
# 675

attr(QC.expr,"thresholds")
#lower   higher 
#546.1979       Inf 

#cells that would be thrown away only for one reason 
sum(QC.mit != QC.expr)
#   1835

colData(sce) <- cbind(colData(sce), df)
sce$discard <- as.logical(QC.mit+QC.expr)

pdf("FIGURES/01_dicard_genes.pdf", width = 4, height = 6)
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

pdf("FIGURES/01_Mito_discard.pdf", width = 8, height = 4)
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

pdf("FIGURES/01_MA_discarded_cells.pdf", width = 6, height = 4)
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[Mt], logFC[Mt], col="dodgerblue", pch=16)
text(abundance[Mt], logFC[Mt], labels = "Mito")
points(abundance["Malat1"], logFC["Malat1"], col="orange", pch=16)
text(abundance["Malat1"], logFC["Malat1"], labels = "Malat1")
dev.off()

MeanAbundance.damaged <- rowMeans(logcounts(sce.damaged))
topGenes.damaged.id <- order(MeanAbundance.damaged,decreasing = TRUE)[1:30]
topGenes.damaged <- rownames(sce.damaged)[topGenes.damaged.id]
#MeanAbundance.damaged[topGenes.damaged.id]
rowData(sce[names(MeanAbundance.damaged[topGenes.damaged.id])])$SYMBOL
# "mt-Atp6" "mt-Co3"  "mt-Co2"  "mt-Co1"  "mt-Cytb" "mt-Nd4"  "mt-Nd2"  "Cst3"    "mt-Nd1"  "Apoe"   
# [11] "Tmsb4x"  "Fth1"    "Actb"    "Tpt1"    "Rpl13"   "mt-Nd5"  "Eef1a1"  "Rplp1"   "Itm2b"   "H2-D1"  
# [21] "Rps8"    "Fau"     "mt-Nd3"  "Rps29"   "H3f3b"   "Ftl1"    "Rpl41"   "Ctsd"    "Rps24"   "Rpl27a" 


sce <- sce[,!sce$discard]
dim(sce)
#  22542 19509

## Save filtered Data to my dir
#saveRDS(sce,"Rdata/sce_filtered.rds")



###
# AFter cleanup


sce2 <- sce[,colData(sce)$Experiment == "2"]
sce1 <- sce[,colData(sce)$Experiment == "1"]


#
hto.counts.s1 <- as.matrix(counts(altExp(sce1)))
hto.counts.s1 <- hto.counts.s1[c("HTO_A0306","HTO_A0307"),]
rownames(hto.counts.s1) <- c("HCC1954","MDA231")

hto.counts.s2 <- as.matrix(counts(altExp(sce2)))
hto.counts.s2 <- hto.counts.s2[c("HTO_A0308","HTO_A0306"),]
rownames(hto.counts.s2) <- c("MDA231","HCC1954")

#logUMI vs logHTO
logUMIs <- log10(colSums(counts(sce1)) +1)
logHTOs <- log10(colSums(hto.counts.s1) +1)

df.plot <- data.frame(logUMIs = logUMIs,
                      logHTOs = logHTOs,
                      compartment = colData(sce1)$SampleName)

pdf("plots/Exp1_HtosVsUmis_afterQC.pdf")
ggplot(df.plot, aes(logUMIs,logHTOs,color = compartment) ) + geom_point() + ggtitle("Exp1, HTOs vs UMI exprs")
dev.off()

logUMIs.2 <- log10(colSums(counts(sce2)) +1)
logHTOs.2 <- log10(colSums(hto.counts.s2) +1)

df.plot.2 <- data.frame(logUMIs.2 = logUMIs.2,
                        logHTOs.2 = logHTOs.2,
                        compartment = colData(sce2)$SampleName)

pdf("plots/Exp2_HtosVsUmis_afterQC.pdf")
ggplot(df.plot.2, aes(logUMIs.2,logHTOs.2,color = compartment) ) + geom_point()  + ggtitle("Exp2, HTOs vs UMI exprs")
dev.off()


#plot HTOs

df.plot.HTOs.1 <- data.frame(log10(t(hto.counts.s1)+1),
                             compartment = colData(sce1)$SampleName)


df.plot.HTOs.2 <- data.frame(log10(t(hto.counts.s2)+1),
                             compartment = colData(sce2)$SampleName)

pdf("plots/HTOs_afterQC.pdf")
ggplot(df.plot.HTOs.1, aes(HCC1954,MDA231,color = compartment) ) + geom_point() + ggtitle("Exp1, HTOs")
ggplot(df.plot.HTOs.2, aes(HCC1954,MDA231,color = compartment) ) + geom_point() + ggtitle("Exp2, HTOs")
dev.off()



#remove those with zeros
hto.counts.s1 <- hto.counts.s1[,!(colSums(hto.counts.s1) == 0)]
hto.counts.s2 <- hto.counts.s2[,!(colSums(hto.counts.s2) == 0)]


dim(hto.counts.s1)
#     2  23680
dim(hto.counts.s2)
#    2   10340


summary(colSums(hto.counts.s1))
#   M Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1      49     645    1220    1469   53148 
summary(colSums(hto.counts.s2))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1      26      66    1380    1003  461243

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

pdf("plots/HTOs_afterQC_confidence.pdf")
ggplot(df.plot.HTOs.1, aes(HCC1954,MDA231,color = Confident) ) + geom_point() + ggtitle("Exp1, HTOs")
ggplot(df.plot.HTOs.2, aes(HCC1954,MDA231,color = Confident) ) + geom_point() + ggtitle("Exp2, HTOs")
dev.off()

#check the assigned ambient profiles...DIFFERENT!
metadata(hash.stats.s1)$ambient
#  HCC1954   MDA231 
#56.93756 33.39975 

metadata(hash.stats.s2)$ambient
#MMDA231  HCC1954 
#17.71685 34.93480 


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


pdf("plots/HTOs_profiles_ambient_estimation.pdf")
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
# 8756 14924 

table(hash.stats.s2$Confident) 
#FALSE  TRUE 
#  5301  5039 

h1.conf <-  hist(log10(hto.counts.s1[1,hash.stats.s1$Confident & hash.stats.s1$Best == 1]), breaks = 30)
h2.conf <-  hist(log10(hto.counts.s1[2,hash.stats.s1$Confident & hash.stats.s1$Best == 2]), breaks = 30)


pdf("plots/HTOs_assignment_exp1.pdf")
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
#20     578    1001    1575    1757   36429 
summary(colSums(hto.confident.s2))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18     213    1043    2759    3318  461243

hash.stats.s2[hash.stats.s2$Total == 461243,]
# Total      Best    Second     LogFC    LogFC2
#<numeric> <integer> <integer> <numeric> <numeric>
#  Stroma_2-CCACCATTCAACGTGT    461243         2        NA   7.28223        NA
#Doublet Confident
#<logical> <logical>
#  Stroma_2-CCACCATTCAACGTGT        NA      TRUE


hto.confident.s2[,"Stroma_2-CCACCATTCAACGTGT"] #apparently cluster
#MDA231 HCC1954 
#1966  459277 

colSums(hto.confident.s1) %>% sort(decreasing = TRUE) %>% head(20)


r.s1 <- rank(-hash.stats.s1$Total)
color_Confident.s1 = c("black","green")[as.factor(hash.stats.s1$Confident )]

r.s2 <- rank(-hash.stats.s2$Total)
color_Confident.s2 = c("black","green")[as.factor(hash.stats.s2$Confident )]

pdf("plots/HashTag_RankExp1.pdf", width = 10, height = 5)
plot(r.s1, hash.stats.s1$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Confidently assigned, Exp1",col=color_Confident.s1 )
dev.off()
pdf("plots/HashTag_RankExp2.pdf", width = 10, height = 5)
plot(r.s2, hash.stats.s2$Total, log="xy", xlab="Rank", ylab="Total HTO count", main="Confidently assigned, Exp2",col=color_Confident.s2 )
dev.off()


#HereExp 1
#total hashcount over 100
#metadata(hash.stats.s1)$ambient %>% sum
hash.stats.s1$Confident[hash.stats.s1$Total < 90] %>% table
#FALSE  TRUE 
#6953   293    ...we do not trust these
hash.stats.s1$Confident[which(hash.stats.s1$Total < 90)] <- FALSE

#####################################################
#total cells Exp1
table(hash.stats.s1$Best[hash.stats.s1$Confident]) 
#  1    2 
#8188 6443

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
hash.stats.s2$Confident[hash.stats.s2$Total < 50] %>% table
#FALSE  TRUE 
#4460   187    ...we do not trust these
hash.stats.s2$Confident[which(hash.stats.s2$Total < 50)] <- FALSE

#also from the HTO rank graph ... the first4 : doublets
#doublets.remove.s2 <- r.s2 %>% sort(decreasing = TRUE) %>% head(5) %>% names
#hash.stats.s2[doublets.remove.s2,]$Confident <- FALSE

#####################################################
#total cells Exp1
table(hash.stats.s2$Best[hash.stats.s2$Confident]) 
#    1    2 
#2015 2837 

rownames(hto.counts.s2)
#1 = "MDA231", 2 =  "HCC1954"

# Hashtag Antibodies
Ht.2 <- c("MDA231", "HCC1954")
#only confident:
hash.stats.s2 <- hash.stats.s2[hash.stats.s2$Confident,]


#IMPORTANT STEP: select only those cells confidently assigned to the wanted conditions
sce2 <- sce2[ , rownames(hash.stats.s2)]
sce2$cell.line <- Ht.2[hash.stats.s2$Best]





sce <- cbind(sce1,sce2) #COmbine back


## Save Data to my dir
#saveRDS(sce,"Rdata/sce_demultiplexed_cleaned.rds")

# Read-in data 
#sce <- readRDS("Rdata/sce_demultiplexed_cleaned.rds")


