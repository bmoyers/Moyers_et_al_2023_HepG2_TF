#!/usr/bin/env R
#Supplemental_4_13_14.R

################################################################################
#This script is used to produce Supplemental Figures 4, 13, and 14
#
#Run under R 4.1.0
#This script takes as arguments:
# outDir: Path where files are to be written
# finalAnnotationsTFs: Supplemental Table 1
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# dataDir: Path to directory containing processed and prepared data from
#     Agarwal et al 2023, provided in folder "Agarwal2023".
# bindingExprData: Model data produced by binding expression scripts. Provided
#     for download as binding_expr_models_results.rds
#
################################################################################



################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(GenomicRanges)
library(ggplot2)
library(memes)
library(universalmotif)
library(Biostrings)




################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


################################################################################
#This function restricts peaks to the central 51BP around a peak's summit, and
#returns a file formatted for easy conversion into genomicRanges.
################################################################################
formatPeaks <- function(this_peaksFile) {
  this_peaks <- read.table(this_peaksFile, header=F, sep="\t", stringsAsFactors=F)
  this_peaks <- as.data.frame(this_peaks)
  this_peaks[,2] <- as.numeric(as.character(this_peaks[,2]))
  this_peaks[,3] <- as.numeric(as.character(this_peaks[,3]))
  this_peaks[,10] <- as.numeric(as.character(this_peaks[,10]))
  newStarts <- this_peaks[,2]+this_peaks[,10]-25
  newEnds <- this_peaks[,2]+this_peaks[,10]+25

  this_peaks[,2] <- newStarts
  this_peaks[,3] <- newEnds
  colnames(this_peaks)[1:3] <- c("chr", "start", "end")
  return(this_peaks)
}


################################################################################
#This function grabs expression levels of named genes.
################################################################################
getNamesAndExpr <- function(thisCounts, thisExpr, theGenes) {

  for (i in 1:nrow(thisExpr)) {
    thisExpr[i,1] <- strsplit(thisExpr[i,1], split=".", fixed=T)[[1]][1]
  }

  thisExpr <- thisExpr[thisExpr[,1]%in%theGenes[,1],]
  theGenes <- theGenes[theGenes[,1]%in%thisExpr[,1],]

  theGenes <- theGenes[theGenes[,3]%in%thisCounts[,1],]
  thisCounts <- thisCounts[thisCounts[,1]%in%theGenes[,3],]

  theLines <- c()
  toRemove <- c()
  for (i in 1:nrow(thisCounts)) {
    thisGene <- thisCounts[i,1]
    thisSet <- theGenes[theGenes[,3]==thisGene,]
    exprSet <- thisExpr[thisExpr[,1]%in%thisSet[,1],]
    if(nrow(exprSet)!=1) {
      toRemove <- append(toRemove, i)
    }
    if(nrow(exprSet)==1) { theLines <- rbind(theLines, c(exprSet))}
  }

  thisCounts <- thisCounts[-toRemove,]

  finalTable <- cbind(thisCounts[,1:4], theLines)

  return(finalTable)
}




################################################################################
################################################################################
#Being Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
finalAnnotationsTFs <- args[2]
exprDir <- args[3]
dataDir <- args[4]
bindingExprData <- args[5]


finalAnnotationsTFs <- read.delim(finalAnnotationsTFs, header=T, sep="\t", stringsAsFactors=F)
finalAnnotationsTFs <- finalAnnotationsTFs[finalAnnotationsTFs[,"Preferred_NonPref"]=="Preferred",]
finalAnnotationsTFs <- finalAnnotationsTFs[finalAnnotationsTFs[,"Final.HA.Annotation"]=="TF",]




################################################################################
################################################################################
#Read in Agarwal's MPRA data.
################################################################################
################################################################################


testTable <- paste(dataDir, "sequence_locations_fullSet.txt", sep="")
testTable <- read.table(testTable, header=T, sep="\t", stringsAsFactors=F)

for (i in 1:nrow(testTable)) {
  if(i%%1000==0) { print(paste(i, nrow(testTable)))}
  testTable[i,11] <- substr(testTable[i,11], 16, 215)
}

testQuants <- paste(dataDir, "sequence_scores_fullSet.txt", sep="")
testQuants <- read.table(testQuants, header=T, sep="\t", stringsAsFactors=F)



################################################################################
################################################################################
#Identify control sequences and signals for easier visualization
################################################################################
################################################################################

#ControlNames
controlNames <- testTable[grep("^C", testTable[,1]),1:2]
controlQuants <- testQuants[testQuants[,1]%in%controlNames[,1],]

controlDF <- c()
for (i in 1:nrow(controlQuants)) {
  thisType <- controlNames[controlNames[,1]==controlQuants[i,1],2]
  controlDF <- rbind(controlDF, c(controlQuants[i,"mean"],thisType))
}
controlDF <- as.data.frame(controlDF)
controlDF[,1] <- as.numeric(as.character(controlDF[,1]))
controlDF[,2] <- as.character(controlDF[,2])
controlDF[,2] <- gsub(", Smith", "", controlDF[,2])
colnames(controlDF) <- c("Scores", "Label")



################################################################################
#Set up the genomicRanges for the mpra data.
################################################################################

testTable <- testTable[testTable$name%in%testQuants$name,]


new_seqDF <- testTable[,c("chr.hg38", "start.hg38", "stop.hg38", "str.hg38", "name")]
colnames(new_seqDF) <- gsub(".hg38", "", colnames(new_seqDF))
colnames(new_seqDF)[4] <- "strand"
new_seqDF <- new_seqDF[!is.na(new_seqDF[,2]),]
new_seqDF_gr <- makeGRangesFromDataFrame(new_seqDF)


################################################################################
################################################################################
#Initial analysis of the number of factors bound.
#determine the number of factors bound within an element region (summit +/-25bp)
#and bin elements by number bound.  Boxplot of expression across those bins.
################################################################################
################################################################################


theExprs <- list.files(exprDir, full.names=T, pattern="Preferred")
theExprs <- theExprs[-grep("H[3|4]K", theExprs)]


boundMatrix <- matrix(nrow=nrow(new_seqDF), ncol=length(theExprs), data=0)
allTFs <- c()
for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisTF <- strsplit(theExprs[i], split="/")[[1]]
  thisTF <- thisTF[length(thisTF)]
  thisTF <- strsplit(thisTF, split="_")[[1]][1]
  thisTF <- gsub("-FLAG", "", thisTF)
  thisTF <- gsub("-eGFP", "", thisTF)
  thesePeaks <- formatPeaks(theExprs[i])
  thesePeaks_gr <- makeGRangesFromDataFrame(thesePeaks)
  theIntersect <- as.data.frame(findOverlaps(thesePeaks_gr, new_seqDF_gr, type="within"))
  boundMatrix[unique(theIntersect[,2]),i] <- 1

  allTFs <- c(allTFs, thisTF)

}
colnames(boundMatrix) <- allTFs
rownames(boundMatrix) <- new_seqDF$name

numBound <- rowSums(boundMatrix)


binNumbers <- c(0, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400)
binNames <- c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300", "301-400", "401+")
theBins <- rep("0", nrow(new_seqDF))
for (i in 1:length(binNumbers)) {
  theBins[numBound>binNumbers[i]] <- binNames[i]
}

graphDF <- cbind(new_seqDF$name, numBound, theBins)
graphDF <- as.data.frame(graphDF)
colnames(graphDF) <- c("name", "preciseNumber", "DAPsBound")
graphDF$numBound <- numBound
graphDF <- graphDF[graphDF$name%in%testQuants[,1],]


theSignal <- rep(NA, nrow(graphDF))
for (i in 1:nrow(graphDF)) {
  theSignal[i] <- testQuants[testQuants[,1]==graphDF[i,"name"],"mean"]
}

graphDF$Signal <- theSignal

graphDF_genes <- graphDF[grep("ENSG", graphDF[,1]),]
graphDF_distal <- graphDF[-grep("ENSG", graphDF[,1]),]


graphDF_genes$category <- rep("Promoter", nrow(graphDF_genes))
graphDF_distal$category <- rep("Distal", nrow(graphDF_distal))

###Add in the controls.

alt_controlDF <- data.frame(name=rep(NA, nrow(controlDF)), preciseNumber=rep(NA, nrow(controlDF)), DAPsBound=controlDF$Label, numBound=rep(NA, nrow(controlDF)), Signal=controlDF$Scores, category=rep("Control", nrow(controlDF)))

this_graphDF <- rbind(graphDF_genes, graphDF_distal, alt_controlDF)
this_graphDF$DAPsBound <- factor(this_graphDF$DAPsBound, levels=c("0", binNames, "negative", "positive"))


category <- c("Control", "Distal", "Promoter")
custom_col <- c("grey", "yellow", "red")


saveFile <- paste(outDir, "Supplemental_4.pdf", sep="")
p <- ggplot(this_graphDF, aes(x=DAPsBound, y=Signal, fill=category)) + theme_bw() + geom_boxplot() + ylab("MPRA Signal") + xlab("DAPs Bound") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.title=element_text(size=25), legend.text=element_text(size=20)) + ylim(-3,5) +
  scale_fill_manual(name="Category", values = custom_col, na.value="grey50")
ggsave(saveFile)




################################################################################
#We now want to look at various metrics of activator and repressor binding.
#We read in the candidates first.  Identify the candidate repressors (should
#be 26 of them), candidate activators (should be 200 of them-- restrict to the
#top 26 based on frac_sig and median_estimate) and some randomly-selected factors
#which fall into neither of those categories (note: set your seed at this point
#for reproducible analyses).
################################################################################


################################################################################
#Read in and identify the candidate repressors.
################################################################################


theData <- readRDS(bindingExprData)

models_df <- as.data.frame(theData[[1]])
summary_df <- as.data.frame(theData[[2]])


candidateDF <- c()
for (i in 2:nrow(summary_df)) {
  print(paste(i, nrow(summary_df)))
  thisTF <- summary_df[i,1]
  this_fraction_significant <- summary_df[i,"frac_sig"]
  #mean_est <- summary_df[i,"mean_est"]
  this_miniModels <- models_df[models_df[,1]==thisTF,c("estimate", "pearson_cor", "spearman_cor", "rsq", "p.value")]
  this_miniModels <- this_miniModels[as.numeric(as.character(this_miniModels[,"p.value"]))<=0.05,]
  median_estimate <- median(as.numeric(as.character(this_miniModels[,"estimate"])))
  #theName <- strsplit(summary_df[i,19], split=" ")[[1]][1]
  theName <- toupper(thisTF)


  miniDF <- cbind(this_miniModels, rep(theName, nrow(this_miniModels)), rep(this_fraction_significant, nrow(this_miniModels)), rep(median_estimate, nrow(this_miniModels)))
  colnames(miniDF) <- c("estimate", "pearson_cor", "spearman_cor", "rsq", "p.value", "TF", "frac_sig", "median_estimate")
  candidateDF <- rbind(candidateDF, miniDF)
}

candidateDF <- as.data.frame(candidateDF)
for (i in c(1,2,3,4,5,7,8)) {candidateDF[,i] <- as.numeric(as.character(candidateDF[,i]))}
candidateDF[,6] <- as.character(candidateDF[,6])

candidateDF[candidateDF[,6]=="NKX31",6] <- "NKX3-1"

candidateRepressors <- candidateDF[,6:8]
candidateRepressors <- unique(candidateRepressors)
candidateRepressors <- candidateRepressors[candidateRepressors[,3]<0,]
candidateRepressors <- candidateRepressors[candidateRepressors[,2]>=0.5,]
candidateRepressors <- candidateRepressors[!is.na(candidateRepressors[,1]),]

Repressors <- candidateRepressors[,1]
length(Repressors[Repressors%in%colnames(boundMatrix)])
#26



candidateActivators <- candidateDF[,6:8]
candidateActivators <- unique(candidateActivators)
candidateActivators <- candidateActivators[candidateActivators[,3]>0,]
candidateActivators <- candidateActivators[candidateActivators[,2]>=0.5,]
candidateActivators <- candidateActivators[!is.na(candidateActivators[,1]),]
candidateActivators <- candidateActivators[order(candidateActivators$frac_sig, candidateActivators$median_estimate, decreasing=T),]

Activators <- candidateActivators[,1]
length(Activators[Activators%in%colnames(boundMatrix)])
#26



candidateNeutrals <- candidateDF[,6:8]
candidateNeutrals <- unique(candidateNeutrals)
candidateNeutrals <- candidateNeutrals[!candidateNeutrals[,1]%in%c(Repressors, Activators),]


Neutrals <- candidateNeutrals[,1]
set.seed(1)
Neutrals <- Neutrals[sample(1:length(Neutrals), 26, replace=F)]

length(Neutrals[Neutrals%in%colnames(boundMatrix)])



################################################################################
#Repressors. For these, specifically, we want to rule out the effect of
#REST in addition to the normal analysis.
################################################################################


boundMatrix_repressors <- boundMatrix[,colnames(boundMatrix)%in%Repressors]
repressorsSums <- rowSums(boundMatrix_repressors)

haveRepressorsBound <- rownames(boundMatrix)[repressorsSums>0]
RepressorBound <- rep("NoRepressorBound", nrow(graphDF))
RepressorBound[graphDF[,1]%in%haveRepressorsBound] <- "RepressorBound"

graphDF$RepressorBound <- RepressorBound

graphDF[as.character(graphDF$DAPsBound)=="negative","RepressorBound"] <- "NegativeControl"
graphDF[as.character(graphDF$DAPsBound)=="positive","RepressorBound"] <- "PositiveControl"

RepressorBound_noRest <- graphDF$RepressorBound
REST_bound <- boundMatrix_repressors <- boundMatrix[,"REST"]
haveRESTBound <- rownames(boundMatrix)[REST_bound>0]
RepressorBound_noRest[which(graphDF$name%in%haveRESTBound)] <- "RESTBound"
graphDF$RepressorBound_noRest <- RepressorBound_noRest

graphDF$RepressorBound <- factor(graphDF$RepressorBound, levels=c("RepressorBound", "NoRepressorBound", "NegativeControl", "PositiveControl"))
graphDF$RepressorBound_noRest <- factor(graphDF$RepressorBound_noRest, levels=c("RepressorBound", "RESTBound", "NoRepressorBound", "NegativeControl", "PositiveControl"))

graphDF_promoters <- graphDF[grep("ENSG", graphDF[,1]),]
graphDF_promoters <- rbind(graphDF_promoters, graphDF[is.na(graphDF[,1]),])



category <- c("RepressorBound", "RestBound", "Unbound", "NegativeControl", "PositiveControl")
alt_custom_col <- c("red", "orange", "purple", "darkgrey", "lightgrey")


saveFile <- paste(outDir, "Supplemental_14.pdf", sep="")
p <- ggplot(graphDF_promoters, aes(x=DAPsBound, y=Signal, fill=RepressorBound_noRest)) + theme_bw() + geom_boxplot() +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20))  +
  ylim(-3,4) + scale_fill_manual(name="Category", values = alt_custom_col, na.value="grey50") + xlab("DAPs Bound") + ylab("MPRA Signal")
ggsave(saveFile)


this_graphDF_promoters <- graphDF_promoters
this_graphDF_distal <- graphDF[-grep("ENSG", graphDF[,1]),]

this_graphDF_promoters_repressors <- this_graphDF_promoters
this_graphDF_distal_repressors <- this_graphDF_distal


################################################################################
#Activators
################################################################################

boundMatrix_activators <- boundMatrix[,colnames(boundMatrix)%in%Activators]
activatorSums <- rowSums(boundMatrix_activators)

haveActivatorsBound <- rownames(boundMatrix)[activatorSums>0]

ActivatorBound <- rep("NoActivatorBound", nrow(graphDF))
ActivatorBound[graphDF[,1]%in%haveActivatorsBound] <- "ActivatorBound"

graphDF$ActivatorBound <- ActivatorBound
graphDF[as.character(graphDF$DAPsBound)=="negative","ActivatorBound"] <- "NegativeControl"
graphDF[as.character(graphDF$DAPsBound)=="positive","ActivatorBound"] <- "PositiveControl"

graphDF$ActivatorBound <- factor(graphDF$ActivatorBound, levels=c("ActivatorBound", "NoActivatorBound", "NegativeControl", "PositiveControl"))

graphDF_promoters <- graphDF[grep("ENSG", graphDF[,1]),]
graphDF_promoters <- rbind(graphDF_promoters, graphDF[is.na(graphDF[,1]),])

graphDF_distal <- graphDF[-grep("ENSG", graphDF[,1]),]


this_graphDF_promoters <- graphDF_promoters
this_graphDF_distal <- graphDF[-grep("ENSG", graphDF[,1]),]

this_graphDF_promoters_activators <- this_graphDF_promoters
this_graphDF_distal_activators <- this_graphDF_distal


################################################################################
#Do the above for each of the individual 26 candidate randomly-selected.
################################################################################

boundMatrix_neutrals <- boundMatrix[,colnames(boundMatrix)%in%Neutrals]
neutralSums <- rowSums(boundMatrix_neutrals)

haveNeutralsBound <- rownames(boundMatrix)[neutralSums>0]

NeutralBound <- rep("NoOtherBound", nrow(graphDF))
NeutralBound[graphDF[,1]%in%haveNeutralsBound] <- "OtherBound"

graphDF$NeutralBound <- NeutralBound
graphDF[as.character(graphDF$DAPsBound)=="negative","NeutralBound"] <- "NegativeControl"
graphDF[as.character(graphDF$DAPsBound)=="positive","NeutralBound"] <- "PositiveControl"

graphDF$NeutralBound <- factor(graphDF$NeutralBound, levels=c("OtherBound", "NoOtherBound", "NegativeControl", "PositiveControl"))

graphDF_promoters <- graphDF[grep("ENSG", graphDF[,1]),]
graphDF_promoters <- rbind(graphDF_promoters, graphDF[is.na(graphDF[,1]),])

graphDF_distal <- graphDF[-grep("ENSG", graphDF[,1]),]

this_graphDF_promoters <- graphDF_promoters
this_graphDF_distal <- graphDF[-grep("ENSG", graphDF[,1]),]

this_graphDF_promoters_random <- this_graphDF_promoters
this_graphDF_distal_random <- this_graphDF_distal


################################################################################
#Make a graph which displays the combination of all of these.
################################################################################

this_repressors <- this_graphDF_promoters_repressors[this_graphDF_promoters_repressors$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300"),c("name", "preciseNumber", "DAPsBound", "Signal", "RepressorBound")]
this_repressors$DAPsBound <- as.character(this_repressors$DAPsBound)
this_repressors <- this_repressors[,1:5]
colnames(this_repressors)[5] <- c("Class")
this_repressors$DAPsBound <- paste("Repressor_", this_repressors$DAPsBound, sep="")
this_repressors$Class <- as.character(this_repressors$Class)
this_repressors[this_repressors$Class=="NoRepressorBound","Class"] <- "Unbound"

this_activators <- this_graphDF_promoters_activators[this_graphDF_promoters_activators$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30"),c("name", "preciseNumber", "DAPsBound", "Signal", "ActivatorBound")]
this_activators$DAPsBound <- as.character(this_activators$DAPsBound)
colnames(this_activators)[5] <- c("Class")
this_activators$DAPsBound <- paste("Activator_", this_activators$DAPsBound, sep="")
this_activators$Class <- as.character(this_activators$Class)
this_activators[this_activators$Class=="NoActivatorBound","Class"] <- "Unbound"

this_random <- this_graphDF_promoters_random[this_graphDF_promoters_random$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200"),c("name", "preciseNumber", "DAPsBound", "Signal", "NeutralBound")]
this_random$DAPsBound <- as.character(this_random$DAPsBound)
colnames(this_random)[5] <- c("Class")
this_random$DAPsBound <- paste("Other_", this_random$DAPsBound, sep="")
this_random$Class <- as.character(this_random$Class)
this_random[this_random$Class=="NoOtherBound","Class"] <- "Unbound"



fun_new_graphDF <- rbind(this_activators, this_repressors, this_random)
fun_new_graphDF$DAPsBound <- factor(fun_new_graphDF$DAPsBound,
  levels=c("Activator_1-5", "Activator_6-10", "Activator_11-20", "Activator_21-30",
  "Repressor_1-5", "Repressor_6-10", "Repressor_11-20", "Repressor_21-30", "Repressor_31-40", "Repressor_41-50", "Repressor_51-100", "Repressor_101-200", "Repressor_201-300",
  "Other_1-5", "Other_6-10", "Other_11-20", "Other_21-30", "Other_31-40", "Other_41-50", "Other_51-100", "Other_101-200"))

#


custom_col=c("lightblue", "purple", "red", "grey")

saveFile <- paste(outDir, "Supplemental_13.pdf", sep="")
p <- ggplot(fun_new_graphDF, aes(x=DAPsBound, y=Signal, fill=Class)) + theme_bw() + geom_boxplot() + ylab("MPRA Signal") + xlab("DAPs Bound") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20))  +
  ylim(-2,3.5) + scale_fill_manual(name="Category", values = custom_col, na.value="grey50") +
  scale_x_discrete(labels=c('1-5', '6-10', '11-20', '21-30', '1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '101-200', '201-300', '1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '101-200'))
ggsave(saveFile, width=10)




#
