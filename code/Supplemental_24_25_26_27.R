#!/usr/bin/env R
#Supplemental_24_25_26_27.R

################################################################################
#This script is used to produce Supplemental Figures 24, 25, 26, and 27
#
#Run under R 3.6.1
#This script takes as arguments:
# outDir: Path where files are to be written
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# exclusionList: Path to the ENCODE exclusion list, under accession ENCFF356LFX
# chromSizes: path to a table of chromosome sizes under hg38, provided as hg38.chrom.sizes
# thePromoters: file containing refseq TSS with gene names +/- 500bp,
#     provided for download as refseq_genes_unique_TSS_1000.bed
# abc: Path to ABC connections, provided by Jessie Engreitz
#
################################################################################

################################################################################
################################################################################
#Load libraries.
################################################################################
################################################################################


library(GenomicRanges)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(Biostrings)


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################

################################################################################
#This function takes a list of bed files and a genomic ranges object of the regions
#we're interested in.  It then creates a matrix of binding for each bed file
#over each region.
################################################################################
createBoundMatrix <- function(allFiles, chromSizes, theSplit="_", windowSize=2000, exclusionList) {
  chromDF <- c()
  for (i in 1:nrow(chromSizes)) {
    thisVec <- seq(1, as.numeric(as.character(chromSizes[i,2])), by=windowSize)
    vec1 <- thisVec[1:(length(thisVec)-1)]
    vec2 <- thisVec[2:length(thisVec)]-1
    thisDF <- cbind(rep(chromSizes[i,1],length(vec1)),vec1,vec2)
    chromDF <- rbind(chromDF, thisDF)
  }
  chromDF <- as.data.frame(chromDF)
  chromDF[,1] <- as.character(chromDF[,1])
  chromDF[,2] <- as.numeric(as.character(chromDF[,2]))
  chromDF[,3] <- as.numeric(as.character(chromDF[,3]))
  colnames(chromDF) <- c("chr", "start", "end")
  theRegions_gr <- makeGRangesFromDataFrame(chromDF, ignore.strand=T, keep.extra.columns=F)

  exclusionList <- as.data.frame(exclusionList)
  colnames(exclusionList)[1:3] <- c("chr", "start", "end")
  exclusionList_gr <- makeGRangesFromDataFrame(exclusionList, ignore.strand=T)


  exclusionIntersect <- as.data.frame(findOverlaps(theRegions_gr, exclusionList_gr))
  if(nrow(exclusionIntersect)>0) {
    chromDF <- chromDF[-unique(exclusionIntersect[,1]),]
    chromDF <- as.data.frame(chromDF)
    colnames(chromDF)[1:3] <- c("chr", "start", "end")
    theRegions_gr <- makeGRangesFromDataFrame(chromDF)
  }

  TFBound <- matrix(nrow=length(theRegions_gr), ncol=length(allFiles), data=0)
  theColnames <- c()
  for (i in 1:length(allFiles)) {
    thisTF <- strsplit(allFiles[i], split="/")[[1]]
    thisTF <- thisTF[length(thisTF)]
    thisTF <- strsplit(thisTF, split=theSplit)[[1]][1]
    thisDF <- read.table(allFiles[i], header=F, sep="\t", stringsAsFactors=F)
    thisDF <- thisDF[thisDF[,1]%in%chromSizes[,1],]
    thisDF <- as.data.frame(thisDF)
    colnames(thisDF)[1:3] <- c("chr", "start", "end")
    thisDF_gr <- makeGRangesFromDataFrame(thisDF)

    thisIntersect <- as.data.frame(findOverlaps(theRegions_gr, thisDF_gr))
    TFBound[unique(thisIntersect[,1]),i] <- 1
    theColnames <- c(theColnames, thisTF)
  }
  colnames(TFBound) <- theColnames
  rownames(TFBound) <- paste(chromDF[,1], ":", chromDF[,2], "-", chromDF[,3], sep="")
  return(TFBound)
}



################################################################################
#This function takes a dataframe used for graphing PCA results
#(first 3 columns should be chr, start, end), a path to a bed file, a factorName
#which will be the new column name, and optional labels.  It then intersects
#the regions in the graphDF with the bedfile to determine which regions
#overlap with the bed of interest, and returns the dataframe with the new column.
################################################################################
add_bed_factor <- function(graphDF, filePath, factorName, nonLabel=NULL, Label=NULL) {
  graphDF <- graphDF[,colnames(graphDF)!=factorName]
  if( (is.null(nonLabel) && !is.null(Label)) || (!is.null(nonLabel) && is.null(Label))) {
    print("Either both labels must be provided, or neither.")
    return(0)
  }
  if(!is.null(dim(filePath))) {
    theFile <- filePath
  }
  if(is.null(dim(filePath))) {
    theFile <- read.delim(filePath, header=F, sep="\t", stringsAsFactors=F)
  }
  thisDF <- as.data.frame(theFile)
  colnames(thisDF)[1:3] <- c("chr", "start", "end")
  thisDF_gr <- makeGRangesFromDataFrame(thisDF)

  graphDF <- as.data.frame(graphDF)
  colnames(graphDF)[1:3] <- c("chr", "start", "end")
  graphDF_gr <- makeGRangesFromDataFrame(graphDF)

  theIntersect <- as.data.frame(findOverlaps(graphDF_gr, thisDF_gr))
  if(is.null(nonLabel)) {
    thisLabel <- rep(paste("Non-", factorName, sep=""), nrow(graphDF))
    thisLabel[unique(theIntersect[,1])] <- factorName
    thisLabel <- factor(thisLabel, levels=c(factorName, paste("Non-", factorName, sep="")))
  }
  if(!is.null(nonLabel)) {
    thisLabel <- rep(nonLabel, nrow(graphDF))
    thisLabel[unique(theIntersect[,1])] <- Label
    thisLabel <- factor(thisLabel, levels=c(Label, nonLabel))
  }
  graphDF <- cbind(graphDF, thisLabel)
  colnames(graphDF)[ncol(graphDF)] <- factorName
  return(graphDF)
}


################################################################################
#This function takes a dataframe used for graphing PCA results
#(first 3 columns should be chr, start, end), and a dataframe in bed format which
#has some numeric factor that we would like to add. It then performs an intersect
#between the two, and assigns a numeric value to each region in the relevant bed file.
################################################################################
add_numeric_factor <- function(graphDF, numericDF, factorName, numericColumn=NA, nonStatus=NA, multipleHandle="max") {
  if(is.na(numericColumn)) {
    print("Column with numeric value must be specified.")
    return(0)
  }
  if(!multipleHandle%in%c("max", "min", "replace")) {
    print("multipleHandle variable must be specified as max, min, or replace.")
    return(0)
  }
  numericDF <- as.data.frame(numericDF)
  colnames(numericDF)[1:3] <- c("chr", "start", "end")
  numericDF_gr <- makeGRangesFromDataFrame(numericDF)

  graphDF <- as.data.frame(graphDF)
  colnames(graphDF)[1:3] <- c("chr", "start", "end")
  graphDF_gr <- makeGRangesFromDataFrame(graphDF)
  theIntersect <- as.data.frame(findOverlaps(graphDF_gr, numericDF_gr))

  thisLabel <- rep(nonStatus, nrow(graphDF))
  for (i in 1:length(theIntersect[,1])) {
    thisVal <- numericDF[theIntersect[i,2],numericColumn]
    if( (is.na(nonStatus) && is.na(thisLabel[theIntersect[i,1]]))  || (!is.na(nonStatus) && thisLabel[theIntersect[i,1]]==nonStatus) ) { thisLabel[theIntersect[i,1]] <- thisVal }
    if( (is.na(nonStatus) && !is.na(thisLabel[theIntersect[i,1]])) || (!is.na(nonStatus) && thisLabel[theIntersect[i,1]]!=nonStatus ) ) {
      if(multipleHandle=="max") { thisLabel[theIntersect[i,1]] <- max(c(thisLabel[theIntersect[i,1]], thisVal)) }
      if(multipleHandle=="min") { thisLabel[theIntersect[i,1]] <- min(c(thisLabel[theIntersect[i,1]], thisVal)) }
      if(multipleHandle=="replace") { thisLabel[theIntersect[i,1]] <- thisVal }
    }
  }
  graphDF <- cbind(graphDF, thisLabel)
  colnames(graphDF)[ncol(graphDF)] <- factorName
  return(graphDF)
}



################################################################################
#This function takes a dataframe with the columns "firstPRC", "secondPRC", and "Measure"
#to plot the two PRCs and color data by measure.  xlabel and ylabel are self-explanatory,
#and the "measure" argument should be a string descriptor of the data plotted.
#In addition to a colored scatterplot, there will also be a density (for factor)
#and box (for continuous) plot above and to the right of the plot showing how
#the measure of interest is distributed across each PC.
#saveFile should be the destination that the file is to be saved to.
################################################################################
graphPC_info <- function(thisDF, saveFile, xlabel, ylabel, measure, alpha=0.5) {

  if(is.factor(thisDF[,"Measure"])) {
    if(1%in%thisDF[,"Measure"]) {
      thisDF$Measure <- as.numeric(as.character(thisDF$Measure))
      thisDF[thisDF[,"Measure"]==1,"Measure"] <- measure
      thisDF[thisDF[,"Measure"]==0,"Measure"] <- paste("Non-", measure, sep="")
      thisDF$Measure <- factor(thisDF$Measure, levels=c(measure, paste("Non-", measure, sep="")))
    }
    mini_thisDF <- thisDF[!is.na(thisDF[,"Measure"]),]
    scatterPlot <- ggplot(mini_thisDF, aes(firstPRC, secondPRC, color=Measure)) +
      xlab(xlabel) + ylab(ylabel) + theme_classic(base_size=15) +
      guides(col = guide_legend(override.aes = list(shape = 15, size = 3))) +
      geom_point() + theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title = element_blank(), legend.text=element_text(size = 5))

    xdensity <- ggplot(mini_thisDF, aes(firstPRC, fill=Measure)) +
      geom_density(alpha=.5) + theme(legend.position = "none") +
      xlab(xlabel) + theme_classic(base_size=10) + theme(legend.position = "none")

    ydensity <- ggplot(mini_thisDF, aes(secondPRC, fill=Measure)) +
      geom_density(alpha=.5) + theme(legend.position = "none") +
      xlab(ylabel) + theme_classic(base_size=10) + coord_flip() + theme(legend.position = "none")
  }


  if(!is.factor(thisDF[,"Measure"])) {
    mini_thisDF <- thisDF[!is.na(thisDF[,"Measure"]),]
    #mini_thisDF <- mini_thisDF[order(mini_thisDF$Measure),]
    breaks1 <- seq(floor(min(mini_thisDF[,1])),ceiling(max(mini_thisDF[,1])), by=1)
    breaks2 <- seq(floor(min(mini_thisDF[,2])),ceiling(max(mini_thisDF[,2])), by=1)
    group1 <- rep(breaks1[1], nrow(mini_thisDF))
    for (l in 2:length(breaks1)) {group1[mini_thisDF[,1]>=breaks1[l]] <- breaks1[l]}
    group2 <- rep(breaks2[1], nrow(mini_thisDF))
    for (l in 2:length(breaks2)) {group2[mini_thisDF[,2]>=breaks2[l]] <- breaks2[l]}
    mini_thisDF$group1 <- group1
    mini_thisDF$group1 <- factor(mini_thisDF$group1, levels=breaks1)
    mini_thisDF$group2 <- group2
    mini_thisDF$group2 <- factor(mini_thisDF$group2, levels=breaks2)

    xMax_1 <- 0
    for (j in 1:(length(breaks1))) {
      thisSet <- mini_thisDF[mini_thisDF$group1==breaks1[j],]
      thisSet <- thisSet$Measure
      thisScore <- quantile(thisSet, 0.90)
      if(!is.na(thisScore)) { xMax_1 <- max(c(xMax_1, thisScore)) }
    }
    xMax_2 <- 0
    for (j in 1:(length(breaks2))) {
      thisSet <- mini_thisDF[mini_thisDF$group2==breaks2[j],]
      thisSet <- thisSet$Measure
      thisScore <- quantile(thisSet, 0.90)
      if(!is.na(thisScore)) { xMax_2 <- max(c(xMax_2, thisScore)) }
    }

    scatterPlot <- ggplot(mini_thisDF, aes(firstPRC, secondPRC, color=Measure)) +
      xlab(xlabel) + ylab(ylabel) + theme_classic(base_size=15) +
      scale_colour_gradient(low = "yellow", high = "black") +
      guides(col = guide_colourbar(barwidth=0.5)) +
      geom_point(alpha=alpha) + theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title = element_blank(), legend.text=element_text(size = 5))

    xdensity <- ggplot(mini_thisDF, aes(x=group1, y=Measure)) +
      geom_boxplot(outlier.shape=NA) + theme(legend.position = "none") + ylim(0, xMax_1) +
      xlab(xlabel) + ylab(measure) + theme_classic(base_size=10) + theme(axis.text.x=element_text(angle = 90)) + theme(legend.position = "none")

    ydensity <- ggplot(mini_thisDF, aes(x=group2, y=Measure)) +
      geom_boxplot(outlier.shape=NA) + theme(legend.position = "none") + ylim(0, xMax_2) +
      xlab(ylabel) + ylab(measure) + theme_classic(base_size=10) + coord_flip() + theme(legend.position = "none")
  }

  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank() ) +
    ggtitle(measure) + theme(plot.title = element_text(face="bold", size=22))
  p <- grid.arrange(xdensity, blankPlot, scatterPlot, ydensity,
    ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

  ggsave(saveFile, plot=p, device = "png", type="cairo")
}



################################################################################
################################################################################
#Begin script.
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
exprDir <- args[2]
exclusionList <- args[3]
chromSizes <- args[4]
thePromoters <- args[5]
abc <- args[6]

################################################################################
#Build the TFBound matrix, and restrict to regions where 2 or more TFs are bound
#for the PCA.
#Run PCA
################################################################################

allFiles <- list.files(exprDir, full.names=T, pattern="Preferred")
allFiles <- allFiles[-grep("H3K[3|4]", allFiles)]




exclusionList <- read.table(exclusionList, header=F, sep="\t", stringsAsFactors=F)
exclusionList <- as.data.frame(exclusionList)
colnames(exclusionList)[1:3] <- c("chr", "start", "end")
exclusionList_gr <- makeGRangesFromDataFrame(exclusionList)



chromSizes <- read.table(chromSizes, header=F, stringsAsFactors=F)
chromSizes <- chromSizes[1:24,]

chromDF <- c()
for (i in 1:nrow(chromSizes)) {
  thisVec <- seq(1, as.numeric(as.character(chromSizes[i,2])), by=2000)
  vec1 <- thisVec[1:(length(thisVec)-1)]
  vec2 <- thisVec[2:length(thisVec)]-1
  thisDF <- cbind(rep(chromSizes[i,1],length(vec1)),vec1,vec2)
  chromDF <- rbind(chromDF, thisDF)
}
chromDF <- as.data.frame(chromDF)
chromDF[,1] <- as.character(chromDF[,1])
chromDF[,2] <- as.numeric(as.character(chromDF[,2]))
chromDF[,3] <- as.numeric(as.character(chromDF[,3]))
colnames(chromDF) <- c("chr", "start", "end")
theRegions_gr <- makeGRangesFromDataFrame(chromDF, ignore.strand=T, keep.extra.columns=F)


exclusionIntersect <- as.data.frame(findOverlaps(theRegions_gr, exclusionList_gr))
if(nrow(exclusionIntersect)>0) {
  chromDF <- chromDF[-unique(exclusionIntersect[,1]),]
  row.names(chromDF) <- c()
  chromDF <- as.data.frame(chromDF)
  colnames(chromDF)[1:3] <- c("chr", "start", "end")
  theRegions_gr <- makeGRangesFromDataFrame(chromDF)
}









TFBound <- createBoundMatrix(allFiles, chromSizes, theSplit="_", windowSize=2000, exclusionList=exclusionList)
colnames(TFBound) <- make.unique(colnames(TFBound), sep="_")

theRowSums <- rowSums(TFBound)
TFBound_reduced <- TFBound[theRowSums>2,]

thisPCA <- princomp(TFBound_reduced)

PoV <- (thisPCA$sdev^2)/sum((thisPCA$sdev)^2)
prc1_lab=paste("PC1 (", round(as.numeric(PoV[1]), digits=4)*100, "% Variance)", sep="")
prc2_lab=paste("PC2 (", round(as.numeric(PoV[2]), digits=4)*100, "% Variance)", sep="")
prc3_lab=paste("PC3 (", round(as.numeric(PoV[3]), digits=4)*100, "% Variance)", sep="")
prc4_lab=paste("PC4 (", round(as.numeric(PoV[4]), digits=4)*100, "% Variance)", sep="")
lab_vec <- c(prc1_lab, prc2_lab, prc3_lab, prc4_lab)



################################################################################
#Construct a dataframe with PCs and various biological measures
#that we would like to explore.
################################################################################

#Set up the initial data frame, including number bound.
theScores <- thisPCA$scores
regions <- matrix(nrow=nrow(TFBound_reduced), ncol=3, data=unlist(strsplit(rownames(TFBound_reduced), split=":|-")), byrow=T)

graphDF <- cbind(regions, theScores[,1:4], rowSums(TFBound_reduced))
graphDF <- as.data.frame(graphDF)
graphDF[,1] <- as.character(graphDF[,1])
for (i in 2:ncol(graphDF)) { graphDF[,i] <- as.numeric(as.character(graphDF[,i]))}
colnames(graphDF) <- c("chr", "start", "end", "prc1", "prc2", "prc3", "prc4", "numBound")
graphDF_gr <- makeGRangesFromDataFrame(graphDF, ignore.strand=T)


#Add promoter status, as measured by whether or not the region overlaps a TSS+/-500bp.
TSS <- thePromoters
graphDF <- add_bed_factor(graphDF, TSS, "Promoter", "Distal", "Promoter")




#Add HepG2 ABC values.

#Add HepG2 ABC values.
abc <- read.delim(abc, header=T, sep="\t", stringsAsFactors=F)


abc <- as.data.frame(abc)
abc <- abc[,c("chr", "start", "end", "ABC.Score")]
abc <- unique(abc)
graphDF <- add_numeric_factor(graphDF, abc, "ABC", numericColumn="ABC.Score", nonStatus=NA, multipleHandle="max")
#rm(abc)




################################################################################
#Quick side-note: make a graph of numBound versus ABC score.
################################################################################

theCor <- cor.test(graphDF$numBound, graphDF$ABC, method="spearman")
theTitle <- paste("ABC v Number Bound, Sp.Rho=", round(theCor$estimate, digits=4), ", p=", round(theCor$p.value, digits=4), sep="")

saveFile <- paste(outDir, "Supplemental_27.pdf", sep="")
p <- ggplot(graphDF, aes(x=ABC, numBound)) +
  theme_classic(base_size=15) + geom_point(alpha=0.1) + ggtitle(theTitle)
ggsave(saveFile)



################################################################################
################################################################################
#Plotting principle components colored by various factors to determine what
#each of the PCs is telling us.
################################################################################
################################################################################


################################################################################
#Plot PC1 and PC2 highlighting each of the qualities.
################################################################################



saveFile <- paste(outDir, "Supplemental_24.png", sep="")
thisDF <- graphDF[,c("prc1", "prc2", "numBound")]
colnames(thisDF) <- c("firstPRC", "secondPRC", "Measure")
graphPC_info(thisDF, saveFile, lab_vec[1], lab_vec[2], "numBound")


saveFile <- paste(outDir, "Supplemental_25.png", sep="")
thisDF <- graphDF[,c("prc1", "prc2", "Promoter")]
colnames(thisDF) <- c("firstPRC", "secondPRC", "Measure")
graphPC_info(thisDF, saveFile, lab_vec[1], lab_vec[2], "Promoter")

saveFile <- paste(outDir, "Supplemental_26.png", sep="")
thisDF <- graphDF[,c("prc1", "prc2", "ABC")]
colnames(thisDF) <- c("firstPRC", "secondPRC", "Measure")
graphPC_info(thisDF, saveFile, lab_vec[1], lab_vec[2], "ABC", alpha=1)






#
