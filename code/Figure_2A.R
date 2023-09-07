#!/usr/bin/env R
#Figure_2A.R

################################################################################
#This script is used to produce Figure 2A.
#
#Run under R 3.6.1
#This script takes as arguments:
# outDir: Path where files are to be written
# theMatrix: Table containing information about which TFs had a peak 1kb upstream
#     of the TSS of a given gene, and the TPM of that gene. Provided for download
#     as Binding_Expression_Matrix.txt
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# abc: Path to ABC connections, provided by Jessie Engreitz
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(ggplot2)
library(matrixStats)
library(MASS)
library(GenomicRanges)


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

makeGraph <- function(theY, theX, saveFile, xLab="", yLab="", mainLab="", theCor=NULL, theCor2=NULL) {
  df <- cbind(as.numeric(theY), as.numeric(theX))
  df <- as.data.frame(df, stringsAsFactors=F)
  df$density <- get_density(df[,1], df[,2], n=30)
  ggplot(df, aes(x=theX, y=theY, color = density)) +
    geom_point() + labs(x=xLab, y=yLab) + theme_classic() +
    theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    geom_smooth(method='lm') + xlim(-3,8) + ylim(-3,9)
  ggsave(saveFile)
  #, title=mainLab, subtitle=paste("rho=", round(theCor, digits=4), ", r=", round(theCor2, digits=4), sep="")
}




fitLM <- function(theTable, saveBase, hist_dataFrame, geneNames, abc) {

  #First, remove any columns in which no TF has bound any promoter.
  theTable <- theTable[,colSums(theTable)>0]

  #Sample the genes to build a model on 70% and test on 30%.
  set.seed(1)
  theSample <- sample(1:nrow(theTable), floor(0.7*nrow(theTable)), replace=F)

  training_Table <- theTable[theSample,]
  test_Table <- theTable[-theSample,]
  test_hist_dataFrame <- hist_dataFrame[-theSample,]
  test_geneNames <- geneNames[-theSample]

  theLM <- lm(TPM ~ . , data=training_Table)

  thePred_lm <- predict(theLM, test_Table[,2:ncol(test_Table)])
  theCor <- cor.test(test_Table[,1], thePred_lm, method="spearman")
  theCor2 <- cor.test(test_Table[,1], thePred_lm, method="pearson")
  print(theCor2)

  saveFile_1 <- paste(saveBase, ".pdf", sep="")
  mainLab <- strsplit(saveBase, split="/")[[1]]
  mainLab <- mainLab[length(mainLab)]
  xLab <- "Observed log(TPM)"
  yLab <- "Predicted log(TPM)"
  makeGraph(test_Table[,1], thePred_lm, saveFile_1, xLab=xLab, yLab=yLab, mainLab=mainLab, theCor=theCor$estimate, theCor2=theCor2$estimate)

  saveFile_2 <- paste(saveBase, ".model.txt", sep="")
  save(theLM, file=saveFile_2)


  a_1 <- which(test_Table$TPM<0)
  a_2 <- which(thePred_lm<0)

  a <- a_1[a_1%in%a_2]
  a_rowSums <- rowSums(test_Table[a,2:ncol(test_Table)])

  print(paste("Number of cases where observed and predicted log(TPM)<0: ", length(a)))


  b_1 <- which(test_Table$TPM<0)
  b_2 <- which(thePred_lm>0)

  b <- b_1[b_1%in%b_2]
  b_rowSums <- rowSums(test_Table[b,2:ncol(test_Table)])
  print(paste("observed log(TPM)<0 and predicted log(TPM)>0 nCases: ", length(b)))
  print(paste("mean number of factors bound at these promoters: ", mean(b_rowSums)))

  d_1 <- which(test_Table$TPM>0)
  d_2 <- which(thePred_lm<0)

  d <- d_1[d_1%in%d_2]

  d_rowSums <- rowSums(test_Table[d,2:ncol(test_Table)])

  print(paste("observed log(TPM)>0 and predicted log(TPM)<0 nCases: ", length(d)))
  print(paste("mean number of factors bound at these promoters: ", mean(d_rowSums)))


  b_colSums <- colSums(test_Table[b,2:ncol(test_Table)])
  b_colFractions <- b_colSums / nrow(test_Table[b,])


  d_colSums <- colSums(test_Table[d,2:ncol(test_Table)])
  d_colFractions <- d_colSums / nrow(test_Table[d,])

  print(rbind(b_colFractions[as.numeric(b_colFractions)>=0.8], d_colFractions[as.numeric(b_colFractions)>=0.8]))


  test_hist_dataFrame <- hist_dataFrame[-theSample,]

  b_histColSums <- colSums(test_hist_dataFrame[b,2:ncol(test_hist_dataFrame)])
  b_histColFractions <- b_histColSums/nrow(test_hist_dataFrame[b,])

  d_histColSums <- colSums(test_hist_dataFrame[d,2:ncol(test_hist_dataFrame)])
  d_histColFractions <- d_histColSums/nrow(test_hist_dataFrame[d,])

  general_histColSums <- colSums(test_hist_dataFrame[,2:ncol(test_hist_dataFrame)])
  general_histColFractions <- general_histColSums/nrow(test_hist_dataFrame)

  print(rbind(b_histColFractions, d_histColFractions, general_histColFractions))

  b_geneNames <- test_geneNames[b]
  print(paste("fraction of obs log(TPM)<0 and pred log(TPM)>0 with abc connection: ", length(b_geneNames[b_geneNames%in%abc$TargetGene])/length(b_geneNames)))


  d_geneNames <- test_geneNames[d]
  print(paste("fraction of obs log(TPM)>0 and pred log(TPM)<0 with abc connection: ", length(d_geneNames[d_geneNames%in%abc$TargetGene])/length(d_geneNames)))

}


################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
theMatrix <- args[2]

theMatrix <- read.delim(theMatrix, header=T, sep="\t", stringsAsFactors=F)
theMatrix <- theMatrix[!is.na(theMatrix[,"TPM"]),]
theMatrix_gr <- makeGRangesFromDataFrame(theMatrix)

geneNames <- theMatrix$name
abc <- read.delim(abc, header=T, sep="\t", stringsAsFactors=F)



exprDir <- "/cluster/home/bmoyers/ENCODE_500Plus/Experiment_Set_2021July6/"
peaks <- list.files(exprDir, full.names=T, pattern="Preferred")
peaks_histones <- peaks[grep("H[3|4]K", peaks)]

hist_dataFrame <- data.frame(name=theMatrix$name)
hist_colnames <- c()
for (i in 1:length(peaks_histones)) {
  this_colname <- strsplit(peaks_histones[i], split="/")[[1]]
  this_colname <- this_colname[length(this_colname)]
  this_colname <- strsplit(this_colname, split="_")[[1]][1]
  hist_colnames <- c(hist_colnames, this_colname)

  thisPeaks <- read.delim(peaks_histones[i], header=F, sep="\t", stringsAsFactors=F)
  colnames(thisPeaks)[1:3] <- c("chr", "start", "end")
  thisPeaks_gr <- makeGRangesFromDataFrame(thisPeaks)

  thisOverlaps <- rep(0, nrow(hist_dataFrame))
  thisIntersect <- as.data.frame(findOverlaps(theMatrix_gr, thisPeaks_gr))
  thisOverlaps[unique(thisIntersect[,1])] <- 1
  hist_dataFrame <- cbind(hist_dataFrame, thisOverlaps)
}
colnames(hist_dataFrame)[2:ncol(hist_dataFrame)] <- hist_colnames



theMatrix <- theMatrix[,5:ncol(theMatrix)]
theMatrix[,1] <- log(theMatrix[,1]+0.1)

saveBase <- paste(outDir, "Moyers_Figure2A", sep="")
fitLM(theMatrix, saveBase, hist_dataFrame, geneNames, abc)


#saveFile_1 <- paste(saveBase, ".pdf", sep="")

#
