#!/usr/bin/env R
#Figure_1B_1C_1D_2C_Supplemental_9_10_16_17_20_21_30_34.R

################################################################################
#This script is used to produce Figures 1B, 1C, 1D, and 2C
#as well as Supplemental Figures 9, 10, 16, 17, 20, 21, 34, and 35.
#
#Run under R 4.1.0
#This program takes the following arguments:
# outDir: Path where files are to be written
# finalAnnotationsTFs: Supplemental Table 1
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# dataDir: Path to directory containing the element representation files,
#     provided for download. Files should begin with the string "Element_representation"
# the_info_table: Supplemental Table 7
# gencode: path to gencode annotations-- gencode.v43.basic.annotation.gtf.gz
# abc: Path to ABC connections, provided by Jessie Engreitz
# thePromoters: file containing refseq TSS with gene names +/- 1kb,
#     provided for download as refseq_genes_unique_TSS.bed
# expr1 and expr2: Path to expression levels of genes in HepG2,
#     accession numbers ENCFF533XPJ and ENCFF321JIT
# merged_Regions: Path to a bed file of all merged ChIP-seq peaks regions,
#     provided for download as All_Merged_Peaks.bed.
# bindingExprData: Model data produced by binding expression scripts. Provided
#     for download as binding_expr_models_results.rds
# nullseq: Path to the "nullseq_generate.py" script from the LS-GKM package.
# nullseq_indices: Path to the nullseq_indices for the hg38 genome build
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
library("biomaRt")

options(scipen=10000)


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################



################################################################################
#presets for a blank theme for plotting.
################################################################################
blank_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"))

#

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
the_info_table <- args[5]
gencode <- args[6]
abc <- args[7]
thePromoters <- args[8]
expr1 <- args[9]
expr2 <- args[10]
merged_Regions <- args[11]
bindingExprData <- args[12]
nullseq <- args[13]
nullseq_indices <- args[14]

finalAnnotationsTFs <- read.delim(finalAnnotationsTFs, header=T, sep="\t", stringsAsFactors=F)
finalAnnotationsTFs <- finalAnnotationsTFs[finalAnnotationsTFs[,"Preferred_NonPref"]=="Preferred",]
finalAnnotationsTFs <- finalAnnotationsTFs[finalAnnotationsTFs[,"Final.HA.Annotation"]=="TF",]


################################################################################
#Read in table with sequence info. Grab control sequence names.
################################################################################


the_info_table <- read.table(the_info_table, header=T, sep="\t", stringsAsFactors=F)
the_info_table[,1] <- gsub(">", "", the_info_table[,1])
the_info_table[,2] <- gsub(">", "", the_info_table[,2])

positiveControls <- the_info_table[grep("positive", the_info_table$fullTitles),"shortTitles"]
positiveControls <- c(paste(positiveControls, "_fwd", sep=""), paste(positiveControls, "_rev", sep=""))

negativeControls <- the_info_table[grep("negative", the_info_table$fullTitles),"shortTitles"]
negativeControls <- c(paste(negativeControls, "_fwd", sep=""), paste(negativeControls, "_rev", sep=""))

################################################################################
################################################################################
#Read in our own MPRA data.
################################################################################
################################################################################

################################################################################
#Read in the data for individual reps.
################################################################################


element_saveFiles <- list.files(dataDir, full.names=T, pattern="Element_representation")
#element_saveFiles <- element_saveFiles[grep("_bothRuns", element_saveFiles)]

element_rep_list <- list()
for (i in 1:length(element_saveFiles)) {
  thisElement <- read.table(element_saveFiles[i], header=T, sep="\t", stringsAsFactors=F)
  element_rep_list[[length(element_rep_list)+1]] <- thisElement
}


################################################################################
#Make barcode distribution plots for each of the replicates.
################################################################################

for (i in 1:3) {
  thisSet <- element_rep_list[[i]]
  saveFile <- paste(outDir, "Supplemental_35_rep", i, ".pdf", sep="")
  p <- ggplot(thisSet, aes(x=n_obs_bc)) + geom_histogram() + theme_classic() + ylab("Frequency") + xlab("Barcodes per element") +
    theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5)) + xlim(0,200)
  ggsave(saveFile)

}


################################################################################
#Make Replicate comparison dotplots for each of the replicate pairs.
################################################################################


for (i in 1:(length(element_rep_list)-1)) {
  firstSet <- element_rep_list[[i]]
  firstSet_negative <- firstSet[firstSet$Element%in%negativeControls,]
  firstSet_positive <- firstSet[firstSet$Element%in%positiveControls,]
  firstSet <- firstSet[!firstSet$Element%in%c(negativeControls, positiveControls),]

  for (j in (i+1):(length(element_rep_list))) {
    secondSet <- element_rep_list[[j]]
    secondSet_negative <- secondSet[secondSet$Element%in%negativeControls,]
    secondSet_positive <- secondSet[secondSet$Element%in%positiveControls,]
    secondSet <- secondSet[!secondSet$Element%in%c(negativeControls, positiveControls),]

    element_mini_repMatrix <- as.data.frame(cbind(firstSet[,"log2"], secondSet[,"log2"]))
    colnames(element_mini_repMatrix) <- c("first", "second")
    negative_mini_repMatrix <- as.data.frame(cbind(firstSet_negative[,"log2"], secondSet_negative[,"log2"]))
    colnames(negative_mini_repMatrix) <- c("first", "second")
    positive_mini_repMatrix <- as.data.frame(cbind(firstSet_positive[,"log2"], secondSet_positive[,"log2"]))
    colnames(positive_mini_repMatrix) <- c("first", "second")
    theCor <- cor.test(element_mini_repMatrix[,1], element_mini_repMatrix[,2], method="pearson")
    theTitle <- paste("Rep", i, " v Rep", j, "; R=", round(theCor$estimate, digits=4), ", p=", round(theCor$p.value, digits=4), sep="")

    saveFile <- paste(outDir, "Supplemental_34_rep", i, "_v_rep", j, "_bothRuns.pdf", sep="")
    p <-  ggplot(element_mini_repMatrix, aes(x=first, y=second)) + theme_classic() + geom_point(alpha=0.3) + geom_point(data=negative_mini_repMatrix, colour = "red") + geom_point(data=positive_mini_repMatrix, colour = "blue") + ggtitle(theTitle) +
      theme(axis.text= element_text(size=9), axis.title=element_text(size=20), axis.text.y=element_text(size=9, angle=90, vjust=0.5, hjust=0.5), axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=0.5)) +
      geom_abline(slope=1, intercept=0, color='purple') + xlab(paste("Rep ", i, sep="")) + ylab(paste("Rep ", j, sep=""))
    ggsave(saveFile)

  }
}


################################################################################
#Summing across reps.
################################################################################


dna_counts <- rowSums(cbind(element_rep_list[[1]][,"dna_counts"], element_rep_list[[2]][,"dna_counts"], element_rep_list[[3]][,"dna_counts"]))
rna_counts <- rowSums(cbind(element_rep_list[[1]][,"rna_counts"], element_rep_list[[2]][,"rna_counts"], element_rep_list[[3]][,"rna_counts"]))

combinedMat <- cbind(element_rep_list[[1]][,"Element"], rna_counts, dna_counts)
colnames(combinedMat) <- c("Element", "rna_counts", "dna_counts")
combinedMat <- as.data.frame(combinedMat)
combinedMat$rna_counts <- as.numeric(as.character(combinedMat$rna_counts))
combinedMat$dna_counts <- as.numeric(as.character(combinedMat$dna_counts))
combinedMat$norm_rna_counts <- combinedMat$rna_counts / (sum(combinedMat$rna_counts)/1000000)
combinedMat$norm_dna_counts <- combinedMat$dna_counts / (sum(combinedMat$dna_counts)/1000000)
combinedMat$ratio <- (combinedMat$norm_rna_counts+0.01) / (combinedMat$norm_dna_counts+0.01)
combinedMat$log2 <- log2(combinedMat$ratio)


colnames(combinedMat)[1] <- "name"

################################################################################
#Make a figure exploring the relationship between activity and GC content
#of the insert.
################################################################################


#Trim off the 15bp adapter sequence on each end of the insert.
sequences <- the_info_table[,"sequence"]
sequences <- substr(sequences, 16, nchar(sequences[1]))
sequences <- substr(sequences, 1, nchar(sequences[1])-15)

#convert to a DNAStringSet
dna_sequences <- DNAStringSet(sequences)
dna_sequences_alphFreq <- alphabetFrequency(dna_sequences)
GC_content <- rowSums(dna_sequences_alphFreq[,c("C", "G")]) / rowSums(dna_sequences_alphFreq)

GC <- rep(NA, nrow(combinedMat))
for (i in 1:nrow(combinedMat)) {
  if(i%%100==0) {
    print(paste(i, nrow(combinedMat)))
  }
  thisGC <- GC_content[which(the_info_table[,"shortTitles"]==gsub("_fwd", "", gsub("_rev", "", combinedMat[i,"name"])))]
  GC[i] <- thisGC
}

combinedMat$GC <- GC

cor.test(combinedMat$log2, combinedMat$GC)



saveFile <- paste(outDir, "Supplemental_X_GC_Activity.pdf", sep="")
p <-  ggplot(combinedMat, aes(x=GC, y=log2)) + theme_classic() + geom_point(alpha=0.3) +
  theme(axis.text= element_text(size=9), axis.title=element_text(size=20), axis.text.y=element_text(size=9, angle=90, vjust=0.5, hjust=0.5), axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=0.5)) +
  xlab("%GC") + ylab("MPRA Activity")
ggsave(saveFile)




################################################################################
#Add full titles.
################################################################################


combinedMat$fullTitles <- the_info_table[match(gsub("_fwd", "", gsub("_rev", "", combinedMat[,1])), the_info_table[,1]),"fullTitles"]


################################################################################
#Identify control sequences and signals for easier visualization
################################################################################

#ControlNames
controls <- combinedMat[grep("ontrol", combinedMat$fullTitles),]
controls_BE <- controls[grep("bindingExpr", controls$fullTitles),]
controls_positive <- controls[grep("positive", controls$fullTitles),]
controls_negative <- controls[grep("negative", controls$fullTitles),]

controls_posNeg <- rbind(controls_positive, controls_negative)
controls_posNeg$label <- rep("NegativeControl", nrow(controls_posNeg))
controls_posNeg[grep("ositive", controls_posNeg$fullTitles),"label"] <- "PositiveControl"

################################################################################
#Remove Controls from combinedMat
################################################################################

combinedMat <- combinedMat[-grep("ontrol", combinedMat$fullTitles),]


################################################################################
#Remove bindingExpression and motifMut from combinedMat, as it is not relevant for
#our following tests.
################################################################################

combinedMat <- combinedMat[-grep("bindingExpr", combinedMat$fullTitles),]
combinedMat <- combinedMat[-grep("motifMut", combinedMat$fullTitles),]


################################################################################
#Set up the genomicRanges for the mpra data.
################################################################################


new_seqDF <- matrix(nrow=nrow(combinedMat), ncol=3, data=NA)
for (i in 1:nrow(combinedMat)) {
  thisSet <- strsplit(combinedMat[i,"fullTitles"], split="_")[[1]][1]
  thisSet <- strsplit(thisSet, split=":|-")[[1]]
  new_seqDF[i,] <- thisSet

}
new_seqDF <- as.data.frame(new_seqDF)
new_seqDF[,2] <- as.numeric(as.character(new_seqDF[,2]))
new_seqDF[,3] <- as.numeric(as.character(new_seqDF[,3]))
colnames(new_seqDF) <- c("chr", "start", "end")
new_seqDF$strand <- rep("+")
new_seqDF[grep("_rev", combinedMat[,1]),"strand"] <- "-"
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
rownames(boundMatrix) <- combinedMat$name

numBound <- rowSums(boundMatrix)


binNumbers <- c(0, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400)
binNames <- c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300", "301-400", "401+")
theBins <- rep("0", nrow(new_seqDF))
for (i in 1:length(binNumbers)) {
  theBins[numBound>binNumbers[i]] <- binNames[i]
}

graphDF <- cbind(combinedMat$name, combinedMat$log2, combinedMat$GC, numBound, theBins)
graphDF <- as.data.frame(graphDF)
colnames(graphDF) <- c("name", "Signal", "GC", "numBound", "DAPsBound")

colnames(graphDF)[2] <- "Signal"


graphDF$numBound <- as.numeric(as.character(graphDF$numBound))
graphDF$Signal <- as.numeric(as.character(graphDF$Signal))
graphDF$GC <- as.numeric(as.character(graphDF$GC))


cor.test(as.numeric(as.character(graphDF$Signal)), as.numeric(as.character(graphDF$numBound)))
cor.test(as.numeric(as.character(graphDF$GC)), as.numeric(as.character(graphDF$numBound)))


saveFile <- paste(outDir, "Supplemental_10_Numbound_Activity.pdf", sep="")
p <-  ggplot(graphDF, aes(x=numBound, y=Signal)) + theme_classic() + geom_point(alpha=0.3) +
  theme(axis.text= element_text(size=9), axis.title=element_text(size=20), axis.text.y=element_text(size=9, angle=90, vjust=0.5, hjust=0.5), axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=0.5)) +
  xlab("FactorsBound") + ylab("MPRA Activity")
ggsave(saveFile)


saveFile <- paste(outDir, "Supplemental_10_Numbound_GC.pdf", sep="")
p <-  ggplot(graphDF, aes(x=GC, y=numBound)) + theme_classic() + geom_point(alpha=0.3) +
  theme(axis.text= element_text(size=9), axis.title=element_text(size=20), axis.text.y=element_text(size=9, angle=90, vjust=0.5, hjust=0.5), axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=0.5)) +
  xlab("GC%") + ylab("FactorsBound")
ggsave(saveFile)



################################################################################
#Identify promoter regions based on gencode annotations.
################################################################################

gencode <- read.table(gencode, header=F, sep="\t", stringsAsFactors=F)
gencode <- gencode[gencode[,3]%in%c("gene", "transcript"),]

gencode_TSS <- c()
for (i in 1:nrow(gencode)) {
  if(gencode[i,7]=="+") {
    thisloc <- gencode[i,4]
  }
  if(gencode[i,7]=="-") {
    thisloc <- gencode[i,5]
  }
  gencode_TSS <- rbind(gencode_TSS, c(gencode[i,1], thisloc, thisloc, gencode[i,7]))
}

gencode_TSS <- as.data.frame(gencode_TSS)
colnames(gencode_TSS) <- c("chr", "start", "end", "strand")
gencode_TSS[,2] <- as.numeric(as.character(gencode_TSS[,2]))-100
gencode_TSS[,3] <- as.numeric(as.character(gencode_TSS[,3]))+100
gencode_TSS <- unique(gencode_TSS)
gencode_TSS_gr <- makeGRangesFromDataFrame(gencode_TSS)

promIntersect <- as.data.frame(findOverlaps(new_seqDF_gr, gencode_TSS_gr))


gencode_TSS_distal <- gencode_TSS
gencode_TSS_distal[,2] <- gencode_TSS_distal[,2]-4900
gencode_TSS_distal[,3] <- gencode_TSS_distal[,3]+4900
gencode_TSS_distal_gr <- makeGRangesFromDataFrame(gencode_TSS_distal)

distalIntersect <- as.data.frame(findOverlaps(new_seqDF_gr, gencode_TSS_distal_gr))


graphDF$Signal <- as.numeric(as.character(graphDF$Signal))
graphDF$numBound <- as.numeric(as.character(graphDF$numBound))

graphDF_genes <- graphDF[unique(promIntersect[,1]),]
graphDF_distal <- graphDF[-unique(distalIntersect[,1]),]
graphDF_distal$DAPsBound <- factor(graphDF_distal$DAPsBound, levels=c(0, binNames))

graphDF_genes$category <- rep("Promoter", nrow(graphDF_genes))
graphDF_distal$category <- rep("Distal", nrow(graphDF_distal))


###Add in the controls.

alt_controlDF <- data.frame(name=controls_posNeg$name, Signal=controls_posNeg$log2, numBound=rep(NA, nrow(controls_posNeg)), GC=rep(NA, nrow(controls_posNeg)), DAPsBound=controls_posNeg$label, category=rep("Control", nrow(controls_posNeg)))
this_graphDF <- rbind(graphDF_genes[], graphDF_distal, alt_controlDF)

this_graphDF$Signal <- as.numeric(as.character(this_graphDF$Signal))
this_graphDF$DAPsBound <- factor(this_graphDF$DAPsBound, levels=c("0", binNames, "NegativeControl", "PositiveControl"))


category <- c("Control", "Distal", "Promoter")
custom_col <- c("grey", "yellow", "red")


###
#We remove cases of the Serpina locus and any case where tiling was used,
#as these run the risk of overrepresentation of certain elements.
###

serpinaNames <- the_info_table[grep("erpin", the_info_table$fullTitles),"shortTitles"]
serpinaNames_fwd_rev <- c(paste(serpinaNames, "_fwd", sep=""), paste(serpinaNames, "_rev", sep=""))

this_graphDF_noSerpina <- this_graphDF[!this_graphDF$name%in%serpinaNames_fwd_rev,]

tilingNames <- the_info_table[grep("iling", the_info_table$fullTitles),"shortTitles"]
tilingNames_fwd_rev <- c(paste(tilingNames, "_fwd", sep=""), paste(tilingNames, "_rev", sep=""))

this_graphDF_notiling <- this_graphDF_noSerpina[!this_graphDF_noSerpina$name%in%tilingNames_fwd_rev,]

saveFile <- paste(outDir, "Moyers_Figure1B.pdf", sep="")
p <- ggplot(this_graphDF_notiling, aes(x=DAPsBound, y=Signal, fill=category)) + theme_classic() + geom_boxplot() + ylab("MPRA Signal") + xlab("DAPs Bound") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.title=element_text(size=25), legend.text=element_text(size=20)) + ylim(-3,5) +
  scale_fill_manual(name="Category", values = custom_col, na.value="grey50")
ggsave(saveFile)

numberTable <- table(this_graphDF_notiling[,c("DAPsBound", "category")])
saveFile <- paste(outDir, "Moyers_Figure1B_table.txt", sep="")
write.table(numberTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numberTable <- as.data.frame(read.table(saveFile, header=T, sep="\t", stringsAsFactors=F))

pvals_ttest <- rep(NA, nrow(numberTable))
pvals_promoter_v_negControl <- rep(NA, nrow(numberTable))
pvals_distal_v_negControl <- rep(NA, nrow(numberTable))
controlSet_tTest <- this_graphDF_notiling[this_graphDF_notiling$DAPsBound=="NegativeControl",]
for (i in 1:nrow(numberTable)) {
  thisSet <- this_graphDF_notiling[this_graphDF_notiling$DAPsBound==rownames(numberTable)[i],]
  thisSet_Promoters <- thisSet[thisSet$category=="Promoter",]
  thisSet_Distal <- thisSet[thisSet$category=="Distal",]
  if(nrow(thisSet_Promoters)>0 && nrow(thisSet_Distal)>0) {
    theTest <- t.test(thisSet_Promoters$Signal, thisSet_Distal$Signal)
    pvals_ttest[i] <- theTest$p.value
  }
  if(nrow(thisSet_Promoters)>0) {
    theTest <- t.test(thisSet_Promoters$Signal, controlSet_tTest$Signal)
    pvals_promoter_v_negControl[i] <- theTest$p.value
  }
  if(nrow(thisSet_Distal)>0) {
    theTest <- t.test(thisSet_Distal$Signal, controlSet_tTest$Signal)
    pvals_distal_v_negControl[i] <- theTest$p.value
  }
}
numberTable$tTestpVal <- pvals_ttest
numberTable$Promoter_v_negative_tTestpVal <- pvals_promoter_v_negControl
numberTable$distal_v_negative_tTestpVal <- pvals_distal_v_negControl

write.table(numberTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)


#
graphDF_promoters <- graphDF[unique(promIntersect[,1]),]
graphDF_promoters <- rbind(graphDF_promoters, alt_controlDF[,1:5])
graphDF_promoters$Signal <- as.numeric(as.character(graphDF_promoters$Signal))
graphDF_promoters$DAPsBound <- factor(graphDF_promoters$DAPsBound, levels=c("0", binNames, "NegativeControl", "PositiveControl"))

graphDF_distal <- graphDF[-unique(distalIntersect[,1]),]
graphDF_distal <- rbind(graphDF_distal, alt_controlDF[,1:5])
graphDF_distal$Signal <- as.numeric(as.character(graphDF_distal$Signal))
graphDF_distal$DAPsBound <- factor(graphDF_distal$DAPsBound, levels=c("0", binNames, "NegativeControl", "PositiveControl"))


################################################################################
################################################################################
#Based on the above, we want to look at the expression of connected genes
#for enhancer sequences.  Do this using ABC data, but then expand upon it
#using Mark's hiccups data (which will require a bit more work).
################################################################################
################################################################################

graphDF_distal <- graphDF[-unique(distalIntersect[,1]),]
new_seqDF_distal <- new_seqDF[-unique(distalIntersect[,1]),]
new_seqDF_distal_gr <- makeGRangesFromDataFrame(new_seqDF_distal)


abc <- read.delim(abc, header=T, sep="\t", stringsAsFactors=F)
#abc <- as.data.frame(abc)

abc <- abc[as.numeric(as.character(abc$distance))>=10000,]

abc[,1] <- as.character(abc[,1])
abc[,2] <- as.numeric(as.character(abc[,2]))
abc[,3] <- as.numeric(as.character(abc[,3]))

abc_gr <- makeGRangesFromDataFrame(abc, ignore.strand=T, keep.extra.columns=F)

theIntersect <- as.data.frame(findOverlaps(new_seqDF_distal_gr, abc_gr))

graphDF_distal_abc <- graphDF_distal[unique(theIntersect[,1]),]
graphDF_distal_noAbc <- graphDF_distal[-unique(theIntersect[,1]),]

graphDF_distal_abc$DAPsBound <- factor(graphDF_distal_abc$DAPsBound, levels=c("0", binNames))
graphDF_distal_noAbc$DAPsBound <- factor(graphDF_distal_noAbc$DAPsBound, levels=c("0", binNames))



graphDF_distal_abc$ABC <- "ABC Supported"
graphDF_distal_noAbc$ABC <- "No ABC Support"


thePromoters <- "/cluster/home/bmoyers/hg38_fastas/refseq_genes_unique_TSS.bed"
promoterBaseFile <- read.delim(thePromoters, header=F, sep="\t", stringsAsFactors=F)
finalCounts1 <- promoterBaseFile[,c(4,5)]
finalCounts1 <- cbind(finalCounts1, promoterBaseFile[,1], rowMeans(promoterBaseFile[,c(2,3)]))
colnames(finalCounts1) <- c("GeneName", "GeneBody", "Chromosome", "TSS")


expr1 <- read.delim(expr1, header=T, sep="\t", stringsAsFactors=F)
expr1 <- expr1[625:nrow(expr1),c(1,3:7)]
expr2 <- read.delim(expr2, header=T, sep="\t", stringsAsFactors=F)
expr2 <- expr2[625:nrow(expr2),c(1,3:7)]


finalExpr <- cbind(expr1[,1:3], c((expr1[,4]+expr2[,4])/2), c((expr1[,5]+expr2[,5])/2), c((expr1[,6]+expr2[,6])/2))
colnames(finalExpr) <- colnames(expr1)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
the_genes <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)

fullCounts1 <- getNamesAndExpr(finalCounts1, finalExpr, the_genes)


################################################################################
#For the ABC promoters and candidate distal elements,
#identify the number bound at each distal element. Bin these as before.
#Identify the gene expression level of each connected promoter.
#Create a boxplot of binned number bound versus expression level.
################################################################################


abc_bound <- rep(0, nrow(abc))
for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisTF <- strsplit(theExprs[i], split="/")[[1]]
  thisTF <- thisTF[length(thisTF)]
  thisTF <- strsplit(thisTF, split="_")[[1]][1]
  thisTF <- gsub("-FLAG", "", thisTF)
  thisTF <- gsub("-eGFP", "", thisTF)
  thesePeaks <- formatPeaks(theExprs[i])
  thesePeaks_gr <- makeGRangesFromDataFrame(thesePeaks)
  theIntersect <- as.data.frame(findOverlaps(thesePeaks_gr, abc_gr, type="within"))

  abc_bound[unique(theIntersect[,2])] <- abc_bound[unique(theIntersect[,2])]+1

}


binNumbers <- c(0, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400)
binNames <- c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300", "301-400", "401+")
abcBins <- rep("0", nrow(abc))
for (i in 1:length(binNumbers)) {
  abcBins[abc_bound>binNumbers[i]] <- binNames[i]
}


abc_genes_graphDF <- c()
for (i in 1:nrow(abc)) {
  thisGene <- abc[i,"TargetGene"]
  if(thisGene%in%fullCounts1[,1]) {
    thisTPM <- as.numeric(fullCounts1[fullCounts1[,1]==thisGene,"TPM"])
    thisSet <- c(as.character(abc[i,"name"]), as.numeric(abc_bound[i]), as.character(abcBins[i]), thisGene, thisTPM)
    abc_genes_graphDF <- rbind(abc_genes_graphDF, thisSet)
  }
}


abc_genes_graphDF <- as.data.frame(abc_genes_graphDF)
colnames(abc_genes_graphDF) <- c("regionName", "numBound", "DAPsBound", "geneName", "TPM")
abc_genes_graphDF$numBound <- as.numeric(as.character(abc_genes_graphDF$numBound))
abc_genes_graphDF$TPM <- as.numeric(abc_genes_graphDF$TPM)
abc_genes_graphDF$logTPM <- log(abc_genes_graphDF$TPM+0.1)
abc_genes_graphDF$DAPsBound <- factor(abc_genes_graphDF$DAPsBound, levels=c("0", binNames))

cor.test(abc_genes_graphDF$numBound, abc_genes_graphDF$logTPM)


saveFile <- paste(outDir, "Moyers_Figure1D.pdf", sep="")
p <- ggplot(abc_genes_graphDF, aes(x=DAPsBound, y=logTPM)) + theme_classic() + geom_boxplot() + xlab("DAPs Bound") + ylab("Connected Gene log(TPM)") + ylim(-2.5,7.5) +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5))
ggsave(saveFile)

countTable <- table(abc_genes_graphDF$DAPsBound)
countTable <- cbind(names(countTable), as.numeric(countTable))
saveFile <- paste(outDir, "Moyers_Figure1D_table.txt", sep="")
write.table(countTable, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
countTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)
countTable <- as.data.frame(countTable)
pvals <- c(NA)
zeroSet <- abc_genes_graphDF[abc_genes_graphDF$DAPsBound==0,]
for (i in 1:length(binNames)) {
  thisSet <- abc_genes_graphDF[abc_genes_graphDF$DAPsBound==binNames[i],]
  thisTest <- t.test(thisSet$logTPM, zeroSet$logTPM, alternative="two.sided")
  pvals <- c(pvals, thisTest$p.value)
}
countTable$pVal <- pvals

write.table(countTable, saveFile, row.names=F, col.names=T, sep="\t", quote=F)



################################################################################
#Separate analysis of regions with numBound versus the number with ABC support.
################################################################################

graphDF_distal$DAPsBound <- factor(graphDF_distal$DAPsBound, levels=c(0, binNames))

quickCheckDF <- data.frame(MPRA_Distal=as.numeric(table(graphDF_distal$DAPsBound)), MPRA_Distal_wABC=as.numeric(table(graphDF_distal_abc$DAPsBound)))
quickCheckDF$FractionRetained <- quickCheckDF$MPRA_Distal_wABC / quickCheckDF$MPRA_Distal
rownames(quickCheckDF) <- c(0, binNames)



allBoundRegions <- read.table(merged_Regions, header=F, sep="\t", stringsAsFactors=F)

binNumbers <- c(0, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400)
binNames <- c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300", "301-400", "401+")
allBoundRegionsBins <- rep("0", nrow(allBoundRegions))
for (i in 1:length(binNumbers)) {
  allBoundRegionsBins[allBoundRegions[,4]>binNumbers[i]] <- binNames[i]
}

allBoundRegions$DAPsBound <- allBoundRegionsBins
allBoundRegions$DAPsBound <- factor(allBoundRegions$DAPsBound, levels=c(0, binNames))

quickCheckDF$allABCConnections_MPRA_or_not <- table(abc_genes_graphDF$DAPsBound)
quickCheckDF$allRegions_MPRA_or_not <- table(allBoundRegions$DAPsBound)


quickCheckDF$other_fraction_retained <- quickCheckDF$allABCConnections_MPRA_or_not / quickCheckDF$allRegions_MPRA_or_not

colnames(quickCheckDF) <- c(colnames(quickCheckDF)[1:3], c("allABCConnections", "allRegions", "all_fraction_retained"))


quickCheckDF[,c("allRegions", "allABCConnections", "all_fraction_retained")]

quickCheckDF$DAPsBound <- rownames(quickCheckDF)
quickCheckDF$DAPsBound <- factor(quickCheckDF$DAPsBound, levels=c(0, binNames))
quickCheckDF[1,"all_fraction_retained"] <- NA

saveFile <- paste(outDir, "Moyers_Figure1C.pdf", sep="")
p <- ggplot(quickCheckDF[quickCheckDF$DAPsBound!=0,], aes(x=DAPsBound, y=all_fraction_retained)) + theme_classic() + geom_bar(stat="identity") + ylab("Fraction with ABC Loop") + xlab("DAPs Bound") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5))
ggsave(saveFile)

################################################################################
#Reviewers would like a control comparison for Figure 1C.
################################################################################



control_mergedRegions <- c()
binNames <- c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300", "301-400", "401+")
for (i in 1:length(binNames)) {
  print(paste(i, length(binNames), binNames[i]))
  thisSet <- allBoundRegions[allBoundRegions$DAPsBound==binNames[i],]
  temp_saveFile <- paste(outDir, "temp_binned_allRegions.txt", sep="")
  write.table(thisSet[,1:3], temp_saveFile, row.names=F, col.names=F, sep="\t", quote=F)
  this_controlBed <- paste(outDir, "Moyers_Figure1C_controlBed_abc_", gsub("-", "_", binNames[i]), ".txt", sep="")

  theCommand <- paste("module load cluster/python/2.7.15; module load cluster/python/2.7-modules; python ", nullseq, " -x 1 -m 1000 -r 1 -o ", this_controlBed, " ", temp_saveFile, " hg38 ", nullseq_indices, sep="")
  system(theCommand)

  this_control <- read.table(this_controlBed, header=F, sep="\t", stringsAsFactors=F)
  this_control <- as.data.frame(this_control)
  colnames(this_control) <- colnames(thisSet[,1:3])
  this_control$binned_numBound <- rep(binNames[i], nrow(this_control))
  control_mergedRegions <- rbind(control_mergedRegions, this_control)
#
}

control_mergedRegions <- as.data.frame(control_mergedRegions)
colnames(control_mergedRegions) <- c("chr", "start", "end", "DAPsBound")
control_mergedRegions_gr <- makeGRangesFromDataFrame(control_mergedRegions)
control_mergedRegions$DAPsBound <- factor(control_mergedRegions$DAPsBound, levels=c(binNames))


control_theIntersect <- as.data.frame(findOverlaps(control_mergedRegions_gr, abc_gr))

control_mergedRegions_abc <- control_mergedRegions[unique(control_theIntersect[,1]),]

control_mergedRegions_abc$DAPsBound <- factor(control_mergedRegions_abc$DAPsBound, levels=c(binNames))

control_mergedRegions_abc$ABC <- "ABC Supported"




control_quickCheckDF <- data.frame(MPRA_Distal=as.numeric(table(control_mergedRegions$DAPsBound)), MPRA_Distal_wABC=as.numeric(table(control_mergedRegions_abc$DAPsBound)))
control_quickCheckDF$FractionRetained <- control_quickCheckDF$MPRA_Distal_wABC / control_quickCheckDF$MPRA_Distal
rownames(control_quickCheckDF) <- c(binNames)


control_quickCheckDF$DAPsBound <- rownames(control_quickCheckDF)
control_quickCheckDF$DAPsBound <- factor(control_quickCheckDF$DAPsBound, levels=c(0, binNames))


saveFile <- paste(outDir, "Supplemental_9.pdf", sep="")
p <- ggplot(control_quickCheckDF, aes(x=DAPsBound, y=FractionRetained)) + theme_classic() + geom_bar(stat="identity") + ylab("Fraction with ABC Loop") + xlab("DAPs Bound") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5)) + ylim(0,0.6)
ggsave(saveFile)



this_control_numTable <- cbind(quickCheckDF[quickCheckDF$DAPsBound!=0,c("allRegions", "allABCConnections", "all_fraction_retained")], control_quickCheckDF[,c("MPRA_Distal", "MPRA_Distal_wABC", "FractionRetained")])
colnames(this_control_numTable) <- c("observed_num", "observed_num_wABC", "observed_fraction", "control_num", "control_num_wABC", "control_fraction")

thePvals_control <- c()
for (i in 1:nrow(this_control_numTable)) {
  thisTable_chiSq <- cbind(rbind((this_control_numTable[i,"observed_num"] - this_control_numTable[i,"observed_num_wABC"]), this_control_numTable[i,"observed_num_wABC"]),
    rbind((this_control_numTable[i,"control_num"] - this_control_numTable[i,"control_num_wABC"]), this_control_numTable[i,"control_num_wABC"]))
  this_chisqTest <- chisq.test(thisTable_chiSq)
  thePvals_control <- c(thePvals_control, this_chisqTest$p.value)
}

this_control_numTable$pvalue <- thePvals_control

saveFile <- paste(outDir, "Moyers_Figure1C_control_table.txt", sep="")
write.table(this_control_numTable, saveFile, row.names=F, col.names=T, sep="\t", quote=F)



################################################################################
#We now want to look at various metrics of activator and repressor binding.
#We read in the candidates first.  Identify the candidate repressors, candidate
#activators (restrict to the top 26 based on frac_sig and median_estimate
#for fair comparison to repressors) and some randomly-selected factors
#which fall into neither of those top categories (note: set your seed at this point
#for reproducible analyses).
################################################################################


################################################################################
#Read in and identify the candidate repressors.
################################################################################


bindingExprData <- readRDS(bindingExprData)

models_df <- as.data.frame(bindingExprData[[1]])
summary_df <- as.data.frame(bindingExprData[[2]])


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
candidateDF[candidateDF[,6]=="CSDA",6] <- "YBX3"

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
Activators <- Activators[1:26]
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
#For each of these 3 subsets, create a plot that separates out factors which are
#bound and which are not bound by the factors of interest and plots
#activity binned by total numBound.
################################################################################



################################################################################
#Repressors. For these, specifically, we want to rule out the effect of
#REST in addition to the normal analysis.
################################################################################

boundMatrix_repressors <- boundMatrix[,colnames(boundMatrix)%in%Repressors]
repressorsSums <- rowSums(boundMatrix_repressors)

this_graphDF <- rbind(graphDF, alt_controlDF[,1:5])

haveRepressorsBound <- rownames(boundMatrix)[repressorsSums>0]
RepressorBound <- rep("NoRepressorBound", nrow(this_graphDF))
RepressorBound[this_graphDF[,1]%in%haveRepressorsBound] <- "RepressorBound"


this_graphDF$RepressorBound <- RepressorBound

this_graphDF[as.character(this_graphDF$DAPsBound)=="NegativeControl","RepressorBound"] <- "NegativeControl"
this_graphDF[as.character(this_graphDF$DAPsBound)=="PositiveControl","RepressorBound"] <- "PositiveControl"

RepressorBound_noRest <- this_graphDF$RepressorBound
REST_bound <- boundMatrix_repressors <- boundMatrix[,"REST"]
haveRESTBound <- rownames(boundMatrix)[REST_bound>0]
RepressorBound_noRest[which(this_graphDF$name%in%haveRESTBound)] <- "RESTBound"
this_graphDF$RepressorBound_noRest <- RepressorBound_noRest

this_graphDF$RepressorBound <- factor(this_graphDF$RepressorBound, levels=c("RepressorBound", "NoRepressorBound", "NegativeControl", "PositiveControl"))
this_graphDF$RepressorBound_noRest <- factor(this_graphDF$RepressorBound_noRest, levels=c("RepressorBound", "RESTBound", "NoRepressorBound", "NegativeControl", "PositiveControl"))


this_graphDF$Signal <- as.numeric(as.character(this_graphDF$Signal))
this_graphDF$DAPsBound <- factor(this_graphDF$DAPsBound, levels=c("0", binNames, "NegativeControl", "PositiveControl"))

this_graphDF_promoters <- this_graphDF[unique(promIntersect[,1]),]
this_graphDF_distal <- this_graphDF[-unique(distalIntersect[,1]),]

this_graphDF_promoters <- rbind(this_graphDF_promoters, this_graphDF[grep("ontrol", this_graphDF$RepressorBound),])
this_graphDF_distal <- rbind(this_graphDF_distal, this_graphDF[grep("ontrol", this_graphDF$RepressorBound),])

this_graphDF_promoters <- this_graphDF_promoters[!this_graphDF_promoters$name%in%serpinaNames_fwd_rev,]
this_graphDF_distal <- this_graphDF_distal[!this_graphDF_distal$name%in%serpinaNames_fwd_rev,]





category <- c("RepressorBound", "RestBound", "Unbound", "NegativeControl", "PositiveControl")
alt_custom_col <- c("red", "orange", "purple", "darkgrey", "lightgrey")

saveFile <- paste(outDir, "Supplemental_17.pdf", sep="")
p_repressors_noRest <- ggplot(this_graphDF_promoters, aes(x=DAPsBound, y=Signal, fill=RepressorBound_noRest)) + theme_classic() + geom_boxplot() +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20))  +
  ylim(-2,4.5) + scale_fill_manual(name="Category", values = alt_custom_col, na.value="grey50") + xlab("DAPs Bound") + ylab("MPRA Signal")
ggsave(saveFile)

numbersTable <- table(this_graphDF_promoters[,c("DAPsBound", "RepressorBound_noRest")])
saveFile <- paste(outDir, "Supplemental_Table_12.txt", sep="")
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numbersTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)


pvals_ttest_repressor_v_rest <- rep(NA, nrow(numbersTable))
pvals_ttest_repressor_v_noRepressor <- rep(NA, nrow(numbersTable))
pvals_ttest_rest_v_noRepressor <- rep(NA, nrow(numbersTable))

for (i in 1:nrow(numbersTable)) {
  thisSet <- this_graphDF_promoters[this_graphDF_promoters$DAPsBound==rownames(numbersTable)[i],]
  thisSet_repressors <- thisSet[thisSet$RepressorBound_noRest=="RepressorBound",]
  thisSet_noRepressors <- thisSet[thisSet$RepressorBound_noRest=="NoRepressorBound",]
  thisSet_rest <- thisSet[thisSet$RepressorBound_noRest=="RESTBound",]

  if(nrow(thisSet_repressors)>1 && nrow(thisSet_rest)>1) {
    theTest <- t.test(thisSet_repressors$Signal, thisSet_rest$Signal)
    pvals_ttest_repressor_v_rest[i] <- theTest$p.value
  }

  if(nrow(thisSet_repressors)>1 && nrow(thisSet_noRepressors)>1) {
    theTest <- t.test(thisSet_repressors$Signal, thisSet_noRepressors$Signal)
    pvals_ttest_repressor_v_noRepressor[i] <- theTest$p.value
  }

  if(nrow(thisSet_rest)>1 && nrow(thisSet_noRepressors)>1) {
    theTest <- t.test(thisSet_rest$Signal, thisSet_noRepressors$Signal)
    pvals_ttest_rest_v_noRepressor[i] <- theTest$p.value
  }

}

numbersTable$pvals_ttest_repressor_v_rest <- pvals_ttest_repressor_v_rest
numbersTable$pvals_ttest_repressor_v_noRepressor <- pvals_ttest_repressor_v_noRepressor
numbersTable$pvals_ttest_rest_v_noRepressor <- pvals_ttest_rest_v_noRepressor
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)

this_graphDF_promoters_repressors <- this_graphDF_promoters
this_graphDF_distal_repressors <- this_graphDF_distal



################################################################################
#Activators
################################################################################

boundMatrix_activators <- boundMatrix[,colnames(boundMatrix)%in%Activators]
activatorSums <- rowSums(boundMatrix_activators)


this_graphDF <- rbind(graphDF, alt_controlDF[,1:5])


haveActivatorsBound <- rownames(boundMatrix)[activatorSums>0]

ActivatorBound <- rep("NoActivatorBound", nrow(this_graphDF))
ActivatorBound[this_graphDF[,1]%in%haveActivatorsBound] <- "ActivatorBound"

this_graphDF$ActivatorBound <- ActivatorBound
this_graphDF[as.character(this_graphDF$DAPsBound)=="NegativeControl","ActivatorBound"] <- "NegativeControl"
this_graphDF[as.character(this_graphDF$DAPsBound)=="PositiveControl","ActivatorBound"] <- "PositiveControl"

this_graphDF$ActivatorBound <- factor(this_graphDF$ActivatorBound, levels=c("ActivatorBound", "NoActivatorBound", "NegativeControl", "PositiveControl"))

this_graphDF$Signal <- as.numeric(as.character(this_graphDF$Signal))
this_graphDF$DAPsBound <- factor(this_graphDF$DAPsBound, levels=c("0", binNames, "NegativeControl", "PositiveControl"))


this_graphDF_promoters <- this_graphDF[unique(promIntersect[,1]),]
this_graphDF_distal <- this_graphDF[-unique(distalIntersect[,1]),]

this_graphDF_promoters <- rbind(this_graphDF_promoters, this_graphDF[grep("ontrol", this_graphDF$ActivatorBound),])
this_graphDF_distal <- rbind(this_graphDF_distal, this_graphDF[grep("ontrol", this_graphDF$ActivatorBound),])

this_graphDF_promoters <- this_graphDF_promoters[!this_graphDF_promoters$name%in%serpinaNames_fwd_rev,]
this_graphDF_distal <- this_graphDF_distal[!this_graphDF_distal$name%in%serpinaNames_fwd_rev,]



this_graphDF_promoters_activators <- this_graphDF_promoters
this_graphDF_distal_activators <- this_graphDF_distal


################################################################################
#Do the above for each of the individual 26 randomly-selected.
################################################################################

boundMatrix_neutrals <- boundMatrix[,colnames(boundMatrix)%in%Neutrals]
neutralSums <- rowSums(boundMatrix_neutrals)

this_graphDF <- rbind(graphDF, alt_controlDF[,1:5])

haveNeutralsBound <- rownames(boundMatrix)[neutralSums>0]

NeutralBound <- rep("NoOtherBound", nrow(this_graphDF))
NeutralBound[this_graphDF[,1]%in%haveNeutralsBound] <- "OtherBound"

this_graphDF$NeutralBound <- NeutralBound
this_graphDF[as.character(this_graphDF$DAPsBound)=="NegativeControl","NeutralBound"] <- "NegativeControl"
this_graphDF[as.character(this_graphDF$DAPsBound)=="PositiveControl","NeutralBound"] <- "PositiveControl"

this_graphDF$NeutralBound <- factor(this_graphDF$NeutralBound, levels=c("OtherBound", "NoOtherBound", "NegativeControl", "PositiveControl"))
colnames(this_graphDF)[6] <- "OtherBound"

this_graphDF$Signal <- as.numeric(as.character(this_graphDF$Signal))
this_graphDF$DAPsBound <- factor(this_graphDF$DAPsBound, levels=c("0", binNames, "NegativeControl", "PositiveControl"))



this_graphDF_promoters <- this_graphDF[unique(promIntersect[,1]),]
this_graphDF_distal <- this_graphDF[-unique(distalIntersect[,1]),]

this_graphDF_promoters <- rbind(this_graphDF_promoters, this_graphDF[grep("ontrol", this_graphDF$NeutralBound),])
this_graphDF_distal <- rbind(this_graphDF_distal, this_graphDF[grep("ontrol", this_graphDF$NeutralBound),])

this_graphDF_promoters <- this_graphDF_promoters[!this_graphDF_promoters$name%in%serpinaNames_fwd_rev,]
this_graphDF_distal <- this_graphDF_distal[!this_graphDF_distal$name%in%serpinaNames_fwd_rev,]



this_graphDF_promoters_random <- this_graphDF_promoters
this_graphDF_distal_random <- this_graphDF_distal





################################################################################
#Combine this information into a single graph.
################################################################################


this_repressors <- this_graphDF_promoters_repressors[this_graphDF_promoters_repressors$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300"),]
this_repressors$DAPsBound <- as.character(this_repressors$DAPsBound)
this_repressors <- this_repressors[,1:6]
colnames(this_repressors)[6] <- c("Class")
this_repressors$DAPsBound <- paste("Repressor_", this_repressors$DAPsBound, sep="")
this_repressors$Class <- as.character(this_repressors$Class)
this_repressors[this_repressors$Class=="NoRepressorBound","Class"] <- "Unbound"

this_activators <- this_graphDF_promoters_activators[this_graphDF_promoters_activators$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100"),]
this_activators$DAPsBound <- as.character(this_activators$DAPsBound)
colnames(this_activators)[6] <- c("Class")
this_activators$DAPsBound <- paste("Activator_", this_activators$DAPsBound, sep="")
this_activators$Class <- as.character(this_activators$Class)
this_activators[this_activators$Class=="NoActivatorBound","Class"] <- "Unbound"

this_random <- this_graphDF_promoters_random[this_graphDF_promoters_random$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200"),]
this_random$DAPsBound <- as.character(this_random$DAPsBound)
colnames(this_random)[6] <- c("Class")
this_random$DAPsBound <- paste("Other_", this_random$DAPsBound, sep="")
this_random$Class <- as.character(this_random$Class)
this_random[this_random$Class=="NoOtherBound","Class"] <- "Unbound"



combined_graphDF <- rbind(this_activators, this_repressors, this_random)
combined_graphDF$DAPsBound <- factor(combined_graphDF$DAPsBound,
  levels=c("Activator_1-5", "Activator_6-10", "Activator_11-20", "Activator_21-30", "Activator_31-40", "Activator_41-50", "Activator_51-100",
  "Repressor_1-5", "Repressor_6-10", "Repressor_11-20", "Repressor_21-30", "Repressor_31-40", "Repressor_41-50", "Repressor_51-100", "Repressor_101-200", "Repressor_201-300",
  "Other_1-5", "Other_6-10", "Other_11-20", "Other_21-30", "Other_31-40", "Other_41-50", "Other_51-100", "Other_101-200"))

#


custom_col=c("lightblue", "purple", "red", "grey")

saveFile <- paste(outDir, "Moyers_Figure2C.pdf", sep="")
p <- ggplot(combined_graphDF, aes(x=DAPsBound, y=Signal, fill=Class)) + theme_classic() + geom_boxplot() + ylab("MPRA Signal") + xlab("DAPs Bound") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20))  +
  ylim(-2,3.5) + scale_fill_manual(name="Category", values = custom_col, na.value="grey50") +
  scale_x_discrete(labels=c('1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '101-200', '201-300', '1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '101-200'))
ggsave(saveFile, width=10)




numbersTable <- table(this_repressors[,c("DAPsBound", "Class")])
saveFile <- paste(outDir, "Moyers_Figure2C_repressorsTable.txt", sep="")
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numbersTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

pvals_ttest <- rep(NA, nrow(numbersTable))

for (i in 1:nrow(numbersTable)) {
  thisSet <- this_repressors[this_repressors$DAPsBound==rownames(numbersTable)[i],]
  thisSet_bound <- thisSet[thisSet$Class=="RepressorBound",]
  thisSet_unbound <- thisSet[thisSet$Class=="Unbound",]

  if(nrow(thisSet_bound)>1 && nrow(thisSet_unbound)>1) {
    theTest <- t.test(thisSet_bound$Signal, thisSet_unbound$Signal)
    pvals_ttest[i] <- theTest$p.value
  }
}

numbersTable$pvals_ttest <- pvals_ttest
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)



numbersTable <- table(this_random[,c("DAPsBound", "Class")])
saveFile <- paste(outDir, "Moyers_Figure2C_randomTable.txt", sep="")
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numbersTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

pvals_ttest <- rep(NA, nrow(numbersTable))
for (i in 1:nrow(numbersTable)) {
  thisSet <- this_random[this_random$DAPsBound==rownames(numbersTable)[i],]
  thisSet_bound <- thisSet[thisSet$Class=="OtherBound",]
  thisSet_unbound <- thisSet[thisSet$Class=="Unbound",]

  if(nrow(thisSet_bound)>1 && nrow(thisSet_unbound)>1) {
    theTest <- t.test(thisSet_bound$Signal, thisSet_unbound$Signal)
    pvals_ttest[i] <- theTest$p.value
  }
}

numbersTable$pvals_ttest <- pvals_ttest
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)


numbersTable <- table(this_activators[,c("DAPsBound", "Class")])
saveFile <- paste(outDir, "Moyers_Figure2C_activatorsTable.txt", sep="")
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numbersTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

pvals_ttest <- rep(NA, nrow(numbersTable))
for (i in 1:nrow(numbersTable)) {
  thisSet <- this_activators[this_activators$DAPsBound==rownames(numbersTable)[i],]
  thisSet_bound <- thisSet[thisSet$Class=="ActivatorBound",]
  thisSet_unbound <- thisSet[thisSet$Class=="Unbound",]

  if(nrow(thisSet_bound)>1 && nrow(thisSet_unbound)>1) {
    theTest <- t.test(thisSet_bound$Signal, thisSet_unbound$Signal)
    pvals_ttest[i] <- theTest$p.value
  }
}

numbersTable$pvals_ttest <- pvals_ttest
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)


################################################################################
#Same kind of combo graph, but this time for distal regions.
################################################################################


this_repressors <- this_graphDF_distal_repressors[this_graphDF_distal_repressors$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200", "201-300"),]
this_repressors$DAPsBound <- as.character(this_repressors$DAPsBound)
this_repressors <- this_repressors[,1:6]
colnames(this_repressors)[6] <- c("Class")
this_repressors$DAPsBound <- paste("Repressor_", this_repressors$DAPsBound, sep="")
this_repressors$Class <- as.character(this_repressors$Class)
this_repressors[this_repressors$Class=="NoRepressorBound","Class"] <- "Unbound"

this_activators <- this_graphDF_distal_activators[this_graphDF_distal_activators$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200"),]
this_activators$DAPsBound <- as.character(this_activators$DAPsBound)
colnames(this_activators)[6] <- c("Class")
this_activators$DAPsBound <- paste("Activator_", this_activators$DAPsBound, sep="")
this_activators$Class <- as.character(this_activators$Class)
this_activators[this_activators$Class=="NoActivatorBound","Class"] <- "Unbound"

this_random <- this_graphDF_distal_random[this_graphDF_distal_random$DAPsBound%in%c("1-5", "6-10", "11-20", "21-30", "31-40", "41-50", "51-100", "101-200"),]
this_random$DAPsBound <- as.character(this_random$DAPsBound)
colnames(this_random)[6] <- c("Class")
this_random$DAPsBound <- paste("Other_", this_random$DAPsBound, sep="")
this_random$Class <- as.character(this_random$Class)
this_random[this_random$Class=="NoOtherBound","Class"] <- "Unbound"



fun_new_graphDF <- rbind(this_activators, this_repressors, this_random)
fun_new_graphDF$DAPsBound <- factor(fun_new_graphDF$DAPsBound,
  levels=c("Activator_1-5", "Activator_6-10", "Activator_11-20", "Activator_21-30", "Activator_31-40", "Activator_41-50", "Activator_51-100", "Activator_101-200",
  "Repressor_1-5", "Repressor_6-10", "Repressor_11-20", "Repressor_21-30", "Repressor_31-40", "Repressor_41-50", "Repressor_51-100", "Repressor_101-200", "Repressor_201-300",
  "Other_1-5", "Other_6-10", "Other_11-20", "Other_21-30", "Other_31-40", "Other_41-50", "Other_51-100", "Other_101-200"))

#


custom_col=c("lightblue", "purple", "red", "grey")

saveFile <- paste(outDir, "Supplemental_16.pdf", sep="")
p <- ggplot(fun_new_graphDF, aes(x=DAPsBound, y=Signal, fill=Class)) + geom_boxplot() + theme_classic() + ylab("MPRA Signal") + xlab("DAPs Bound") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20))  +
  ylim(-2,3.5) + scale_fill_manual(name="Category", values = custom_col, na.value="grey50") +
  scale_x_discrete(labels=c('1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '101-200', '1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '101-200', '201-300', '1-5', '6-10', '11-20', '21-30', '31-40', '41-50', '51-100', '101-200'))
ggsave(saveFile, width=10)



numbersTable <- table(this_repressors[,c("DAPsBound", "Class")])
saveFile <- paste(outDir, "Supplemental_Table_11_repressors.txt", sep="")
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numbersTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

pvals_ttest <- rep(NA, nrow(numbersTable))

for (i in 1:nrow(numbersTable)) {
  thisSet <- this_repressors[this_repressors$DAPsBound==rownames(numbersTable)[i],]
  thisSet_bound <- thisSet[thisSet$Class=="RepressorBound",]
  thisSet_unbound <- thisSet[thisSet$Class=="Unbound",]

  if(nrow(thisSet_bound)>1 && nrow(thisSet_unbound)>1) {
    theTest <- t.test(thisSet_bound$Signal, thisSet_unbound$Signal)
    pvals_ttest[i] <- theTest$p.value
  }
}

numbersTable$pvals_ttest <- pvals_ttest
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)



numbersTable <- table(this_random[,c("DAPsBound", "Class")])
saveFile <- paste(outDir, "Supplemental_Table_11_other.txt", sep="")
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numbersTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

pvals_ttest <- rep(NA, nrow(numbersTable))
for (i in 1:nrow(numbersTable)) {
  thisSet <- this_random[this_random$DAPsBound==rownames(numbersTable)[i],]
  thisSet_bound <- thisSet[thisSet$Class=="OtherBound",]
  thisSet_unbound <- thisSet[thisSet$Class=="Unbound",]

  if(nrow(thisSet_bound)>1 && nrow(thisSet_unbound)>1) {
    theTest <- t.test(thisSet_bound$Signal, thisSet_unbound$Signal)
    pvals_ttest[i] <- theTest$p.value
  }
}

numbersTable$pvals_ttest <- pvals_ttest
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)


numbersTable <- table(this_activators[,c("DAPsBound", "Class")])
saveFile <- paste(outDir, "Supplemental_Table_11_activators.txt", sep="")
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
numbersTable <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

pvals_ttest <- rep(NA, nrow(numbersTable))
for (i in 1:nrow(numbersTable)) {
  thisSet <- this_activators[this_activators$DAPsBound==rownames(numbersTable)[i],]
  thisSet_bound <- thisSet[thisSet$Class=="ActivatorBound",]
  thisSet_unbound <- thisSet[thisSet$Class=="Unbound",]

  if(nrow(thisSet_bound)>1 && nrow(thisSet_unbound)>1) {
    theTest <- t.test(thisSet_bound$Signal, thisSet_unbound$Signal)
    pvals_ttest[i] <- theTest$p.value
  }
}

numbersTable$pvals_ttest <- pvals_ttest
write.table(numbersTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)



################################################################################
################################################################################
#For each of the activators and repressors... Along with a subset of other elements,
#namely MAFF, MAFFK, ATF7, JUND, CREB1, FOXA2, HLF, ONECUT1, CEBPA, CEBPB,
#HNF4A, ZNF770, HNRNPK, ZNF121, CEBPG, RBM22, NFIL3, and ZBTB40,
#do the following:
#1) Identify all locations at which the factor is bound.
#2) Identify all locations at which the factor is not bound.
#3) For all of the cases of n total elements bound, match the bound element
#   with an unbound element of the same number of factors bound, with as
#   close as possible to the same identity.
#4) Store the absolute values of the element's activity level, as well as
#   the control element's activity level, and their difference. (obs-cont).
#
#In the end, we would like to make a large boxplot of these elements.
################################################################################
################################################################################

extraElements <- c("MAFF", "MAFK", "ATF7", "JUND", "CREB1", "FOXA2", "HLF", "ONECUT1",
  "CEBPA", "CEBPB", "HNF4A", "ZNF770", "HNRNPK", "ZNF121", "CEBPG", "RBM22", "NFIL3", "ZBTB40")

all_of_interest <- c(Activators, extraElements, Repressors)

boundMatrix_promoters <- boundMatrix[unique(promIntersect[,1]),]
boundMatrix_promoters <- boundMatrix_promoters[!rownames(boundMatrix_promoters)%in%serpinaNames_fwd_rev,]
#which(colnames(boundMatrix_promoters)=="CSDA")
colnames(boundMatrix_promoters)[which(colnames(boundMatrix_promoters)=="CSDA")] <- "YBX3"

all_paired_observations <- c()
the_ttest_table <- c()
for (l in 1:length(all_of_interest)) {

  label <- "Other"
  if(all_of_interest[l]%in%Activators) { label <- "Activator"}
  if(all_of_interest[l]%in%Repressors) { label <- "Repressor"}

  boundMatrix_factor <- boundMatrix_promoters[boundMatrix_promoters[,all_of_interest[l]]==1,]
  boundMatrix_noFactor <- boundMatrix_promoters[boundMatrix_promoters[,all_of_interest[l]]==0,]

  factorSums <- rowSums(boundMatrix_factor)
  noFactorSums <- rowSums(boundMatrix_noFactor)

  uniqueNumbers <- unique(factorSums)
  uniqueNumbers <- uniqueNumbers[order(uniqueNumbers)]

  set.seed(1)
  for (i in uniqueNumbers) {
    print(paste(i, max(uniqueNumbers), all_of_interest[l], l, length(all_of_interest)))

    this_boundSet <- rbind(boundMatrix_factor[factorSums==i,])
    rownames(this_boundSet) <- rownames(boundMatrix_factor)[which(factorSums==i)]
    this_unboundSet <- rbind(boundMatrix_noFactor[noFactorSums==i,])
    rownames(this_unboundSet) <- rownames(boundMatrix_noFactor)[which(noFactorSums==i)]

    original_bound_nrow <- nrow(this_boundSet)
    original_unbound_nrow <- nrow(this_unboundSet)

    if(i==1) {
      allSamples <- this_unboundSet[sample(1:nrow(this_unboundSet), nrow(this_boundSet), replace=F),]
      j <- nrow(this_boundSet)


      this_controlSet <- combinedMat[combinedMat$name%in%rownames(allSamples),]
      this_observedSet <- combinedMat[combinedMat$name%in%rownames(this_boundSet[1:j,]),]

      thisDF <- data.frame(TF=rep(all_of_interest[l], nrow(this_observedSet)), label=rep(label, nrow(this_observedSet)), obsSignal=this_observedSet$log2, compSignal=this_controlSet$log2, obsGC=this_observedSet$GC, compGC=this_controlSet$GC)
      thisDF$difference <- thisDF$obsSignal - thisDF$compSignal
      all_paired_observations <- rbind(all_paired_observations, thisDF)
    }

    if(i>1 && nrow(this_unboundSet)>0 && nrow(this_boundSet)>0) {
      allSamples <- c()

      if(original_bound_nrow<original_unbound_nrow) {

        for (j in 1:nrow(this_boundSet)) {
          if(nrow(rbind(this_unboundSet))>0) {
            difMat <- this_unboundSet
            for (k in 1:ncol(this_boundSet)) {
              difMat[,k] <- difMat[,k]-this_boundSet[j,k]
            }
            difMat <- difMat*difMat
            allDifs <- rowSums(difMat)
            sampleOptions <- which(allDifs==min(allDifs))
            thisSample <- sample(names(sampleOptions), 1)
            allSamples <- rbind(allSamples, this_unboundSet[thisSample,])
            rownames(allSamples)[j] <- thisSample
            this_unboundSet <- this_unboundSet[rownames(this_unboundSet)!=thisSample,]
          }
        }

        this_controlSet <- combinedMat[combinedMat$name%in%rownames(allSamples),]
        this_observedSet <- combinedMat[combinedMat$name%in%rownames(this_boundSet)[1:j],]

      } else {

        for (j in 1:nrow(this_unboundSet)) {
          if(nrow(rbind(this_boundSet))>0) {
            difMat <- rbind(this_boundSet)
            for (k in 1:ncol(this_unboundSet)) {
              difMat[,k] <- difMat[,k]-this_unboundSet[j,k]
            }
            difMat <- difMat*difMat
            allDifs <- rowSums(difMat)
            sampleOptions <- which(allDifs==min(allDifs))
            thisSample <- sample(names(sampleOptions), 1)
            allSamples <- rbind(allSamples, this_boundSet[thisSample,])
            rownames(allSamples)[j] <- thisSample
            theRownames <- rownames(this_boundSet)
            this_boundSet <- rbind(this_boundSet[rownames(this_boundSet)!=thisSample,])
            rownames(this_boundSet) <- theRownames[!theRownames%in%thisSample]
          }
        }

        this_controlSet <- combinedMat[combinedMat$name%in%rownames(this_unboundSet)[1:j],]
        this_observedSet <- combinedMat[combinedMat$name%in%rownames(allSamples),]

      }

      thisDF <- data.frame(TF=rep(all_of_interest[l], nrow(this_observedSet)), label=rep(label, nrow(this_observedSet)), obsSignal=this_observedSet$log2, compSignal=this_controlSet$log2, obsGC=this_observedSet$GC, compGC=this_controlSet$GC)
      thisDF$difference <- thisDF$obsSignal - thisDF$compSignal
      all_paired_observations <- rbind(all_paired_observations, thisDF)
    }

  }

  thisSet <- all_paired_observations[all_paired_observations[,1]==all_of_interest[l],]
  this_ttest <- t.test(thisSet$obsSignal, thisSet$compSignal, paired=T)
  the_ttest_table <- rbind(the_ttest_table, c(all_of_interest[l], label, nrow(thisSet), as.numeric(this_ttest$estimate), as.numeric(this_ttest$p.value)))

}



the_ttest_table <- as.data.frame(the_ttest_table)
colnames(the_ttest_table) <- c("TF", "PredictedType", "nObs", "estimate", "pval")
the_ttest_table$estimate <- as.numeric(the_ttest_table$estimate)
the_ttest_table$pval <- as.numeric(the_ttest_table$pval)
the_ttest_table$nObs <- as.numeric(the_ttest_table$nObs)

the_ttest_table$label <- rep("", nrow(the_ttest_table))
the_ttest_table[the_ttest_table$pval<=0.05,"label"] <- "*"
the_ttest_table[the_ttest_table$pval<=0.0001,"label"] <- "**"
the_ttest_table[the_ttest_table$pval<=2.2e-16,"label"] <- "***"


####Need to reorder these within each category to be from highest to lowest median.
####Also need to reorder as Activator, Repressor, Other.

new_all_of_interest <- c()

act_set <- all_paired_observations[all_paired_observations[,2]=="Activator",]
unique_act_tfs <- unique(act_set[,1])
act_df <- c()
for (i in 1:length(unique_act_tfs)) {
  thisVal <- median(act_set[act_set[,1]==unique_act_tfs[i],"difference"], rm.na=TRUE)
  act_df <- rbind(act_df, c(unique_act_tfs[i], thisVal))
}
act_df <- as.data.frame(act_df)
act_df[,2] <- as.numeric(as.character(act_df[,2]))
act_df <- act_df[order(act_df[,2], decreasing=T),]
new_all_of_interest <- c(new_all_of_interest, act_df[,1])


rep_set <- all_paired_observations[all_paired_observations[,2]=="Repressor",]
unique_rep_tfs <- unique(rep_set[,1])
rep_df <- c()
for (i in 1:length(unique_rep_tfs)) {
  thisVal <- median(rep_set[rep_set[,1]==unique_rep_tfs[i],"difference"], rm.na=TRUE)
  rep_df <- rbind(rep_df, c(unique_rep_tfs[i], thisVal))
}
rep_df <- as.data.frame(rep_df)
rep_df[,2] <- as.numeric(as.character(rep_df[,2]))
rep_df <- rep_df[order(rep_df[,2], decreasing=T),]
new_all_of_interest <- c(new_all_of_interest, rep_df[,1])


oth_set <- all_paired_observations[all_paired_observations[,2]=="Other",]
unique_oth_tfs <- unique(oth_set[,1])
oth_df <- c()
for (i in 1:length(unique_oth_tfs)) {
  thisVal <- median(oth_set[oth_set[,1]==unique_oth_tfs[i],"difference"], rm.na=TRUE)
  oth_df <- rbind(oth_df, c(unique_oth_tfs[i], thisVal))
}
oth_df <- as.data.frame(oth_df)
oth_df[,2] <- as.numeric(as.character(oth_df[,2]))
oth_df <- oth_df[order(oth_df[,2], decreasing=T),]
new_all_of_interest <- c(new_all_of_interest, oth_df[,1])


all_paired_observations$TF <- factor(all_paired_observations$TF, levels=new_all_of_interest)

all_paired_observations$GC_Difference <- all_paired_observations$obsGC - all_paired_observations$compGC

the_ttest_table$color <- rep("red", nrow(the_ttest_table))

against_exp <- c("ZNF687", "ZBTB12", "TBX2", "SALL1", "MGA", "ZFP91", "ZNF501", "ZFX", "ZFY", "ASH2L")
the_ttest_table[the_ttest_table$TF%in%against_exp,"color"] <- "#71716F"

custom_col <- c("lightblue", "purple", "red")

saveFile <- paste(outDir, "Supplemental_20.pdf", sep="")
p <- ggplot(all_paired_observations, aes(x=TF, y=difference, fill=label)) + theme_classic() + geom_boxplot(outlier.shape=NA) + ylim(-2,2) +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
  scale_fill_manual(name="Category", values = custom_col, na.value="grey50") + ylab("Paired Signal Difference") +
  geom_hline(yintercept=0) + geom_text(data=data.frame(), aes(x=the_ttest_table$TF, y=rep(-2, nrow(the_ttest_table)), label=the_ttest_table$label), col=the_ttest_table$color, size=10, angle = 90) #+ coord_flip()
ggsave(saveFile, width=18)

saveFile <- paste(outDir, "Supplemental_20_allData.txt", sep="")
write.table(all_paired_observations, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#all_paired_observations <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

saveFile <- paste(outDir, "Supplemental_Table_15.txt", sep="")
write.table(the_ttest_table, saveFile, row.names=F, col.names=T, sep="\t", quote=F)







saveFile <- paste(outDir, "Supplemental_21_GC_content_difference.pdf", sep="")
p <- ggplot(all_paired_observations, aes(x=GC_Difference, fill=label)) + theme_classic() + geom_density()  + facet_wrap(vars(label)) +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20), strip.text = element_text(size = 20))
ggsave(saveFile, width=18)


#saveFile <- paste(outDir, "Supplemental_X_GC_content_difference.pdf", sep="")
#p <- ggplot(all_paired_observations, aes(x=GC_Difference, fill=label)) + blank_theme + geom_density()  + facet_wrap(vars(label)) +
#  theme(strip.text = element_text(size = 20))
#ggsave(saveFile, width=18)


summary(all_paired_observations[all_paired_observations$label=="Activator","GC_Difference"])
summary(all_paired_observations[all_paired_observations$label=="Repressor","GC_Difference"])
summary(all_paired_observations[all_paired_observations$label=="Other","GC_Difference"])








#




#
