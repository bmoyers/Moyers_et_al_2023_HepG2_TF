#!/usr/bin/env R
#Supplemental_5.R


################################################################################
#This script is used to produce Supplemental Figure 5
#
#Run under R 4.1.0
#This script takes as arguments:
# outDir: Path where files are to be written
# exclusionList: Path to the ENCODE exclusion list, under accession ENCFF356LFX
# cCREs: Path to the Version 4 cCREs provided by Jill Moore in June 2022.
# ATAC_seq: Path to ATAC-seq peaks, under accession ENCFF439EIO
# merged_Regions: Path to a bed file of all merged ChIP-seq peaks regions,
#     provided for download as All_Merged_Peaks.bed.
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
#
################################################################################



################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################


library(GenomicRanges)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(viridis)
require(scales)

################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


blank_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"))

#






################################################################################
################################################################################
#Begin script.
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
exclusionList <- args[2]
cCREs <- args[3]
ATAC_seq <- args[4]
merged_Regions <- args[5]
exprDir <- args[6]

exclusionList <- read.table(exclusionList, header=F, sep="\t", stringsAsFactors=F)
colnames(exclusionList)[1:3] <- c("chr", "start", "end")
exclusionList_gr <- makeGRangesFromDataFrame(exclusionList, ignore.strand=T, keep.extra.columns=F)


cCREs <- read.table(cCREs, header=F, sep="\t", stringsAsFactors=F)
cCREs <- cCREs[,c(1:3,6)]
colnames(cCREs)[1:4] <- c("chr", "start", "end", "state")
cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=TRUE, keep.extra.columns=F)
theStates <- unique(cCREs[,4])

ATAC_seq <- read.delim(ATAC_seq, header=F, sep="\t", stringsAsFactors=F)
ATAC_seq <- ATAC_seq[,1:3]
colnames(ATAC_seq)[1:3] <- c("chr", "start", "end")
ATAC_seq_gr <- makeGRangesFromDataFrame(ATAC_seq, ignore.strand=TRUE, keep.extra.columns=F)


merged_Regions <- read.table(merged_Regions, header=F, sep="\t", stringsAsFactors=F)
colnames(merged_Regions)[1:3] <- c("chr", "start", "end")
merged_Regions_gr <- makeGRangesFromDataFrame(merged_Regions, ignore.strand=T, keep.extra.columns=F)


merged_Intersect <- as.data.frame(findOverlaps(merged_Regions_gr, exclusionList_gr))
if(nrow(merged_Intersect)>0) {
  merged_Regions <- merged_Regions[-unique(merged_Intersect[,1]),]
  merged_Regions_gr <- makeGRangesFromDataFrame(merged_Regions, ignore.strand=T, keep.extra.columns=F)
}


cCREs_Intersect <- as.data.frame(findOverlaps(cCREs_gr, exclusionList_gr))
if(nrow(cCREs_Intersect)>0) {
  cCREs <- cCREs[-unique(cCREs_Intersect[,1]),]
  cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=T, keep.extra.columns=F)
}


ATAC_seq_Intersect <- as.data.frame(findOverlaps(ATAC_seq_gr, exclusionList_gr))
if(nrow(ATAC_seq_Intersect)>0) {
  ATAC_seq <- ATAC_seq[-unique(ATAC_seq_Intersect[,1]),]
  ATAC_seq_gr <- makeGRangesFromDataFrame(ATAC_seq, ignore.strand=T, keep.extra.columns=F)
}





theStates <- c("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "TF", "CA", "None")
custom_col <- c("#F51313", "#F59913", "#F5C813", "#EEA4A0", "#43B6F3", "#56D59C", "#6EA38B", "#808080", "#000000")


thisIntersect <- as.data.frame(findOverlaps(merged_Regions_gr, cCREs_gr))
thisIntersect <- thisIntersect[order(thisIntersect[,1]),]


theExprs <- list.files(exprDir, full.names=T, pattern="Preferred")
theExprs <- theExprs[-grep("H3K[3|4]", theExprs)]



################################################################################
#Here, we're going to create a bound matrix of cCREs inside of
#ATAC-seq peaks.  Subsample various amounts, and determine the
#information gained per TF.
#We'll start with a graph of this information, but depending onthe results
#We may calculate the asymptote.
################################################################################


ATAC_Intersect <- as.data.frame(findOverlaps(cCREs_gr, ATAC_seq_gr))
if(nrow(ATAC_Intersect)>0) {
  cCREs <- cCREs[unique(ATAC_Intersect[,1]),]
  cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=T, keep.extra.columns=F)
}


ATAC_Intersect <- as.data.frame(findOverlaps(merged_Regions_gr, ATAC_seq_gr))
if(nrow(ATAC_Intersect)>0) {
  merged_Regions <- merged_Regions[unique(ATAC_Intersect[,1]),]
  merged_Regions_gr <- makeGRangesFromDataFrame(merged_Regions, ignore.strand=T, keep.extra.columns=F)
}



theExprs <- list.files(exprDir, full.names=T, pattern="Preferred")
theExprs <- theExprs[-grep("H3K[3|4]", theExprs)]

boundMat <- matrix(nrow=nrow(cCREs), ncol=length(theExprs), data=0)

for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisFile <- read.delim(theExprs[i], header=F, sep="\t", stringsAsFactors=F)
  colnames(thisFile)[1:3] <- c("chr", "start", "end")
  thisFile_gr <- makeGRangesFromDataFrame(thisFile, ignore.strand=T)

  theIntersect <- as.data.frame(findOverlaps(cCREs_gr, thisFile_gr))

  boundMat[unique(theIntersect[,1]),i] <- 1
}


graphDF_unbound <- cbind(rep(0, 20), rep(0,20))
theSDs <- c(rep(0,20))
for (i in seq(1, ncol(boundMat))) {
  print(paste(i, ncol(boundMat)))
  thePercs <- c()
  for (j in 1:20) {
    thisSample <- sample(1:ncol(boundMat), i, replace=F)
    miniSet <- cbind(boundMat[,thisSample])
    theSums <- rowSums(miniSet)
    unbound <- length(theSums[theSums==0])/length(theSums)
    thePercs <- c(thePercs, unbound)
  }
  thisSD <- sd(thePercs)
  theSDs <- c(theSDs, thisSD)
  graphDF_unbound <- rbind(graphDF_unbound, cbind(thePercs, rep(i, 20)))

}

graphDF_unbound <- as.data.frame(graphDF_unbound)
colnames(graphDF_unbound) <- c("PercentUnbound", "NumDAPs")
graphDF_unbound[,1] <- as.numeric(as.character(graphDF_unbound[,1]))
graphDF_unbound[,2] <- factor(graphDF_unbound[,2], levels=seq(0, max(graphDF_unbound[,2])))
#graphDF_unbound[,2] <- as.numeric(as.character(graphDF_unbound[,2]))



new_graphDF_unbound <- c()
for (i in 1:max(as.numeric(as.character(graphDF_unbound[,2])))) {
  thisSet <- graphDF_unbound[as.numeric(as.character(graphDF_unbound[,2]))==i,]
  thisSet[,2] <- as.numeric(as.character(thisSet[,1]))
  thisMean <- mean(thisSet[,1])
  thisMin <- min(thisSet[,1])
  thisMax <- max(thisSet[,1])
  new_graphDF_unbound <- rbind(new_graphDF_unbound, c(i, thisMean, thisMin, thisMax))
}

new_graphDF_unbound <- as.data.frame(new_graphDF_unbound)
new_graphDF_unbound[,1] <- as.numeric(as.character(new_graphDF_unbound[,1]))
new_graphDF_unbound[,2] <- as.numeric(as.character(new_graphDF_unbound[,2]))
new_graphDF_unbound[,3] <- as.numeric(as.character(new_graphDF_unbound[,3]))
new_graphDF_unbound[,4] <- as.numeric(as.character(new_graphDF_unbound[,4]))

colnames(new_graphDF_unbound) <- c("NumDAPs", "PercentUnbound", "UnboundMin", "UnboundMax")





################################################################################
#Fitting a GLM to this data to predict out to 1096 experiments.
################################################################################

graphDF_unbound_noZero <- as.data.frame(graphDF_unbound[graphDF_unbound$NumDAPs!=0,])

graphDF_unbound_noZero$InverseNumDAPs <- 1/graphDF_unbound_noZero$NumDAPs
graphDF_unbound_noZero$InverseNumDAPsSquared <- (1/graphDF_unbound_noZero$NumDAPs)*(1/graphDF_unbound_noZero$NumDAPs)
graphDF_unbound_noZero$InverseSqrtNumDAPs <- 1/sqrt(graphDF_unbound_noZero$NumDAPs)

intercept <- 1

a3 <- lm(PercentUnbound~., graphDF_unbound_noZero)

predictDF <- data.frame(NumDAPs=0:1639)
predictDF$InverseNumDAPs <- 1/predictDF$NumDAPs
predictDF$InverseNumDAPsSquared <- (1/predictDF$NumDAPs)*(1/predictDF$NumDAPs)
predictDF$InverseSqrtNumDAPs <- 1/sqrt(predictDF$NumDAPs)


b3 <- predict(a3, predictDF)


modelDF <- data.frame(PercentUnbound=b3[2:length(b3)], NumDAPs=1:1639)


saveFile <- paste(outDir, "Supplemental_5.pdf", sep="")
p <- ggplot(new_graphDF_unbound, aes(x=NumDAPs, y=PercentUnbound)) + geom_point() +
  geom_errorbar(aes(ymin=UnboundMin, ymax=UnboundMax), width=.1) + blank_theme +
  geom_line(data=modelDF, aes(x=NumDAPs, y=PercentUnbound), color='red') + geom_vline(xintercept=1096, color="blue") +
  theme(axis.text.x=element_text(size=20, angle=90), axis.text.y=element_text(size=20), axis.title=element_text(size=25)) +
  xlab("Number of DAPs") + ylab("Fraction of cCREs Unbound") + labs(title="") + xlim(0,1200)
ggsave(saveFile)




#
