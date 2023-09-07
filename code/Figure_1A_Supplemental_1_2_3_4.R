#!/usr/bin/env R
#Figure_1A_Supplemental_1_2_3_4.R

################################################################################
#This script is used to produce Figure 1A as well as Supplemental Figures 1, 2, 3, and 4
#
#Run under R 3.6.1
#This script takes as arguments:
# outDir: Path where files are to be written
# exclusionList: Path to the ENCODE exclusion list, under accession ENCFF356LFX
# cCREs: Path to the Version 4 cCREs provided by Jill Moore in June 2022.
# ATAC_seq: Path to ATAC-seq peaks, under accession ENCFF439EIO
# merged_Regions: Path to a bed file of all merged ChIP-seq peaks regions,
#     provided for download as All_Merged_Peaks.bed.
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# nullseq: Path to the "nullseq_generate.py" script from the LS-GKM package.
# nullseq_indices: Path to the nullseq_indices for the hg38 genome build
#
################################################################################

################################################################################
################################################################################
#Load Libraries
#Set Options
################################################################################
################################################################################

options(scipen=10000)

library(GenomicRanges)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(viridis)
require(scales)



################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
exclusionList <- args[2]
cCREs <- args[3]
ATAC_seq <- args[4]
merged_Regions <- args[5]
exprDir <- args[6]
nullseq <- args[7]
nullseq_indices <- args[8]

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


ATAC_Intersect <- as.data.frame(findOverlaps(cCREs_gr, ATAC_seq_gr))
if(nrow(ATAC_Intersect)>0) {
  cCREs <- cCREs[unique(ATAC_Intersect[,1]),]
  cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=T, keep.extra.columns=F)
}





theExprs <- list.files(exprDir, full.names=T, pattern="Preferred")
theExprs <- theExprs[-grep("H[3|4]K", theExprs)]

numBound <- rep(0, nrow(cCREs))

for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisFile <- read.delim(theExprs[i], header=F, sep="\t", stringsAsFactors=F)
  colnames(thisFile)[1:3] <- c("chr", "start", "end")
  thisFile_gr <- makeGRangesFromDataFrame(thisFile, ignore.strand=T)

  theIntersect <- as.data.frame(findOverlaps(cCREs_gr, thisFile_gr))

  numBound[unique(theIntersect[,1])] <- numBound[unique(theIntersect[,1])]+1
}


cCREs$numBound <- numBound


theStates <- c("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "TF", "CA")
custom_col <- c("#F51313", "#F59913", "#F5C813", "#EEA4A0", "#43B6F3", "#56D59C", "#6EA38B", "#808080")



binned_numBound <- c()
for (i in 1:nrow(cCREs)) {
  if(cCREs[i,"numBound"]==0) { binned_numBound <- c(binned_numBound, "0") }
  if(cCREs[i,"numBound"]>0 && cCREs[i,"numBound"]<= 10 ) { binned_numBound <- c(binned_numBound, "1-10") }
  if(cCREs[i,"numBound"]>10 && cCREs[i,"numBound"]<= 25 ) { binned_numBound <- c(binned_numBound, "11-25") }
  if(cCREs[i,"numBound"]>25 && cCREs[i,"numBound"]<= 50 ) { binned_numBound <- c(binned_numBound, "26-50") }
  if(cCREs[i,"numBound"]>50 && cCREs[i,"numBound"]<= 100 ) { binned_numBound <- c(binned_numBound, "51-100") }
  if(cCREs[i,"numBound"]>100 && cCREs[i,"numBound"]<= 200 ) { binned_numBound <- c(binned_numBound, "101-200") }
  if(cCREs[i,"numBound"]>200 && cCREs[i,"numBound"]<= 300 ) { binned_numBound <- c(binned_numBound, "201-300") }
  if(cCREs[i,"numBound"]>300 && cCREs[i,"numBound"]<= 400 ) { binned_numBound <- c(binned_numBound, "301-400") }
  if(cCREs[i,"numBound"]>400 && cCREs[i,"numBound"]<= 500 ) { binned_numBound <- c(binned_numBound, "401-500") }
  if(cCREs[i,"numBound"]>500 && cCREs[i,"numBound"]<= 600 ) { binned_numBound <- c(binned_numBound, "501-600") }
  if(cCREs[i,"numBound"]>600 ) { binned_numBound <- c(binned_numBound, "601+") }

}

cCREs$binned_numBound <- binned_numBound

cCREs$binned_numBound  <- factor(binned_numBound, levels=c("0", "1-10", "11-25", "26-50", "51-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601+"))

graphDF <- c()
for (i in c("0", "1-10", "11-25", "26-50", "51-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601+")) {
  thisSet <- cCREs[cCREs$binned_numBound==i,]
  typeNumbers <- c()
  for (j in 1:length(theStates)) {
    typeNumbers <- c(typeNumbers, nrow(thisSet[thisSet$state==theStates[j],]))
  }

  if(sum(typeNumbers)>0) { typeFractions <- typeNumbers/sum(typeNumbers) }


  graphDF <- rbind(graphDF, cbind(rep(i, length(theStates)+1),
    c(nrow(thisSet), typeNumbers), c(1, typeFractions), c("numSites", theStates)))

}

graphDF <- as.data.frame(graphDF)
graphDF[,1] <- factor(graphDF[,1], levels=c("0", "1-10", "11-25", "26-50", "51-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601+"))
graphDF[,2] <- as.numeric(as.character(graphDF[,2]))
graphDF[,3] <- as.numeric(as.character(graphDF[,3]))
graphDF[,4] <- as.character(graphDF[,4])
graphDF[,4] <- factor(graphDF[,4], levels=theStates)
colnames(graphDF) <- c("NumBound", "NumSites", "FractionSites", "Type")
graphDF <- graphDF[!is.na(graphDF$Type),]

saveFile <- paste(outDir, "Supplemental_3.pdf", sep="")
ggplot(data=graphDF[graphDF$Type!="numSites",], aes(x=NumBound, y=NumSites, fill=Type)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("DAPS Bound") + ylab("Number of cCREs") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)

saveFile <- paste(outDir, "Supplemental_4.pdf", sep="")
ggplot(data=graphDF[graphDF$Type!="numSites",], aes(x=NumBound, y=FractionSites, fill=Type)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("DAPs Bound") + ylab("Fraction of cCREs") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)

################################################################################
#Need to make Matched control sequences for the above.
#For each bin of numBound, get matched sequence. Then perform overlap with
#the cCRE regions to determine what number of each nullseq falls into the
#cCRE categories.
################################################################################

control_cCREs <- c()
bins <- c("0", "1-10", "11-25", "26-50", "51-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601+")
for (i in 1:length(bins)) {
  print(paste(i, length(bins), bins[i]))
  thisSet <- cCREs[cCREs$binned_numBound==bins[i],]
  temp_saveFile <- paste(outDir, "temp_binned_cCREs.txt", sep="")
  write.table(thisSet[,1:3], temp_saveFile, row.names=F, col.names=F, sep="\t", quote=F)
  this_controlBed <- paste(outDir, "Supplemental_3_controlBed_binned_", gsub("-", "_", bins[i]), ".txt", sep="")

  theCommand <- paste("module load cluster/python/2.7.15; module load cluster/python/2.7-modules; python ", nullseq, " -x 1 -m 1000 -r 1 -o ", this_controlBed, " ", temp_saveFile, " hg38 ", nullseq_indices, sep="")
  system(theCommand)

  this_control <- read.table(this_controlBed, header=F, sep="\t", stringsAsFactors=F)
  this_control <- as.data.frame(this_control)
  colnames(this_control) <- colnames(thisSet[,1:3])
  this_control$binned_numBound <- rep(bins[i], nrow(this_control))
  control_cCREs <- rbind(control_cCREs, this_control)
#
}

#Need to overlap with cCREs, because the goal is to compare these control regions'
#cCRE distribution with the observed distribution.

control_cCREs_gr <- makeGRangesFromDataFrame(control_cCREs)
theIntersect <- as.data.frame(findOverlaps(control_cCREs_gr, cCREs_gr))


new_control_cCREs <- c()
for (i in 1:nrow(theIntersect)) {
  if(i%%100==0) {print(paste(i, nrow(theIntersect)))}
  thisLine <- control_cCREs[theIntersect[i,1],]
  thisType <- cCREs[theIntersect[i,2],"state"]
  thisLine$state <- thisType
  new_control_cCREs <- rbind(new_control_cCREs, thisLine)
}



alt_graphDF <- c()
for (i in c("0", "1-10", "11-25", "26-50", "51-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601+")) {
  thisSet <- new_control_cCREs[new_control_cCREs$binned_numBound==i,]
  typeNumbers <- c()
  for (j in 1:length(theStates)) {
    typeNumbers <- c(typeNumbers, nrow(thisSet[thisSet$state==theStates[j],]))
  }

  if(sum(typeNumbers)>0) { typeFractions <- typeNumbers/sum(typeNumbers) }


  alt_graphDF <- rbind(alt_graphDF, cbind(rep(i, length(theStates)+1),
    c(nrow(thisSet), typeNumbers), c(1, typeFractions), c("numSites", theStates)))

}

alt_graphDF <- as.data.frame(alt_graphDF)
alt_graphDF[,1] <- factor(alt_graphDF[,1], levels=c("0", "1-10", "11-25", "26-50", "51-100", "101-200", "201-300", "301-400", "401-500", "501-600", "601+"))
alt_graphDF[,2] <- as.numeric(as.character(alt_graphDF[,2]))
alt_graphDF[,3] <- as.numeric(as.character(alt_graphDF[,3]))
alt_graphDF[,4] <- as.character(alt_graphDF[,4])
alt_graphDF[,4] <- factor(alt_graphDF[,4], levels=theStates)
colnames(alt_graphDF) <- c("NumBound", "NumSites", "FractionSites", "Type")
alt_graphDF <- alt_graphDF[!is.na(alt_graphDF$Type),]


saveFile <- paste(outDir, "Supplemental_3_controlComparison.pdf", sep="")
ggplot(data=alt_graphDF[alt_graphDF$Type!="numSites",], aes(x=NumBound, y=NumSites, fill=Type)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("DAPS Bound") + ylab("Number of cCREs") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)


saveFile <- paste(outDir, "Supplemental_4_controlComparison.pdf", sep="")
ggplot(data=alt_graphDF[alt_graphDF$Type!="numSites",], aes(x=NumBound, y=FractionSites, fill=Type)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("DAPs Bound") + ylab("Fraction of cCREs") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)


pval_table <- c()
for (i in 1:length(bins)) {
  thisSet_observed <- graphDF[graphDF$NumBound==bins[i],]
  thisSet_control <- alt_graphDF[alt_graphDF$NumBound==bins[i],]
  this_chi_table <- cbind(thisSet_observed$NumSites, thisSet_control$NumSites)
  this_chi_table <- this_chi_table[rowSums(this_chi_table)>0,]
  theTest <- chisq.test(this_chi_table)
  pval_table <- rbind(pval_table, c(bins[i], theTest$p.value))
}
pval_table <- as.data.frame(pval_table)
colnames(pval_table) <- c("bin", "pvalue")
pval_table$pvalue <- as.numeric(as.character(pval_table$pvalue))

saveFile <- paste(outDir, "Supplemental_3_table.txt", sep="")
write.table(pval_table, saveFile, row.names=F, col.names=T, sep="\t", quote=F)


################################################################################
#Now Figure 1A, proportion of each cCRE which is bound and not.
################################################################################


graphDF <- c()
for (i in theStates) {
  thisSet <- cCREs[cCREs$state==i,]
  this_bound <- nrow(thisSet[thisSet$numBound>0,])
  this_unbound <- nrow(thisSet[thisSet$numBound==0,])
  graphDF <- rbind(graphDF, c(i, this_bound, this_bound/(this_bound+this_unbound), "Bound"), c(i, this_unbound, this_unbound/(this_bound+this_unbound), "Unbound"))


}

graphDF <- graphDF[graphDF[,1]!="None",]
graphDF <- as.data.frame(graphDF)
graphDF[,1] <- factor(graphDF[,1], levels=theStates)
graphDF[,2] <- as.numeric(as.character(graphDF[,2]))
graphDF[,3] <- as.numeric(as.character(graphDF[,3]))
graphDF[,4] <- factor(graphDF[,4], levels=c("Unbound", "Bound"))
colnames(graphDF) <- c("State", "Number", "Fraction", "Status")



theStatus <- c("Unbound", "Bound")
#custom_col <- c("purple", "blue")
custom_col <- c("white", "grey")



saveFile <- paste(outDir, "Moyers_Figure1A.pdf", sep="")
ggplot(data=graphDF, aes(x=State, y=Number, fill=Status)) + geom_bar(stat="identity", color="black") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Number of cCREs") + scale_fill_manual(name="", values = custom_col, na.value="grey50") + labs(title="")
ggsave(saveFile)



saveFile <- paste(outDir, "Supplemental_1.pdf", sep="")
ggplot(data=graphDF, aes(x=State, y=Fraction, fill=Status)) + geom_bar(stat="identity", color="black") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Fraction of cCREs") + scale_fill_manual(name="", values = custom_col, na.value="grey50") + labs(title="")
ggsave(saveFile)






################################################################################
#Need to make Matched control sequences for the above.
#For each cCRE set, get matched control sequences, and then determine
#the number and fraction of each category's control which is bound or unbound.
################################################################################

control_cCREs <- c()

for (i in 1:length(theStates)) {
  print(paste(i, length(theStates), theStates[i]))
  thisSet <- cCREs[cCREs$state==theStates[i],]
  temp_saveFile <- paste(outDir, "temp_binned_cCREs.txt", sep="")
  write.table(thisSet[,1:3], temp_saveFile, row.names=F, col.names=F, sep="\t", quote=F)
  this_controlBed <- paste(outDir, "Supplemental_4_controlBed_binned_", theStates[i], ".txt", sep="")

  theCommand <- paste("module load cluster/python/2.7.15; module load cluster/python/2.7-modules; python ", nullseq, " -x 1 -m 1000 -r 1 -o ", this_controlBed, " ", temp_saveFile, " hg38 ", nullseq_indices, sep="")
  system(theCommand)

  this_control <- read.table(this_controlBed, header=F, sep="\t", stringsAsFactors=F)
  this_control <- as.data.frame(this_control)
  colnames(this_control) <- colnames(thisSet[,1:3])
  this_control$state <- rep(theStates[i], nrow(this_control))
  control_cCREs <- rbind(control_cCREs, this_control)
#
}

#Need to determine, for each of these control regions,
#the number of DAPs bound.

control_cCREs_gr <- makeGRangesFromDataFrame(control_cCREs)

theExprs <- list.files(exprDir, full.names=T, pattern="Preferred")
theExprs <- theExprs[-grep("H[3|4]K", theExprs)]

numBound_control <- rep(0, nrow(control_cCREs))

for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisFile <- read.delim(theExprs[i], header=F, sep="\t", stringsAsFactors=F)
  colnames(thisFile)[1:3] <- c("chr", "start", "end")
  thisFile_gr <- makeGRangesFromDataFrame(thisFile, ignore.strand=T)

  theIntersect <- as.data.frame(findOverlaps(control_cCREs_gr, thisFile_gr))

  numBound_control[unique(theIntersect[,1])] <- numBound_control[unique(theIntersect[,1])]+1
}

control_cCREs$numBound <- numBound_control






alt_graphDF <- c()
for (i in theStates) {
  thisSet <- control_cCREs[control_cCREs$state==i,]
  this_bound <- nrow(thisSet[thisSet$numBound>0,])
  this_unbound <- nrow(thisSet[thisSet$numBound==0,])
  alt_graphDF <- rbind(alt_graphDF, c(i, this_bound, this_bound/(this_bound+this_unbound), "Bound"), c(i, this_unbound, this_unbound/(this_bound+this_unbound), "Unbound"))
}

alt_graphDF <- alt_graphDF[alt_graphDF[,1]!="None",]
alt_graphDF <- as.data.frame(alt_graphDF)
alt_graphDF[,1] <- factor(alt_graphDF[,1], levels=theStates)
alt_graphDF[,2] <- as.numeric(as.character(alt_graphDF[,2]))
alt_graphDF[,3] <- as.numeric(as.character(alt_graphDF[,3]))
alt_graphDF[,4] <- factor(alt_graphDF[,4], levels=c("Unbound", "Bound"))
colnames(alt_graphDF) <- c("State", "Number", "Fraction", "Status")



theStatus <- c("Unbound", "Bound")
custom_col <- c("white", "grey")



saveFile <- paste(outDir, "Moyers_Figure1A_controlComparison.pdf", sep="")
ggplot(data=alt_graphDF, aes(x=State, y=Number, fill=Status)) + geom_bar(stat="identity", color="black") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Number of cCREs") + scale_fill_manual(name="", values = custom_col, na.value="grey50") + labs(title="")
ggsave(saveFile)

saveFile <- paste(outDir, "Supplemental_2.pdf", sep="")
ggplot(data=alt_graphDF, aes(x=State, y=Fraction, fill=Status)) + geom_bar(stat="identity", color="black") + theme_classic() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.title=element_text(size=25), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Fraction of cCREs") + scale_fill_manual(name="", values = custom_col, na.value="grey50") + labs(title="")
ggsave(saveFile)



pval_table <- c()
for (i in 1:length(theStates)) {
  thisSet_observed <- graphDF[graphDF$State==theStates[i],]
  thisSet_control <- alt_graphDF[alt_graphDF$State==theStates[i],]
  this_chi_table <- cbind(thisSet_observed$Number, thisSet_control$Number)
  this_chi_table <- this_chi_table[rowSums(this_chi_table)>0,]
  theTest <- chisq.test(this_chi_table)
  pval_table <- rbind(pval_table, c(theStates[i], theTest$p.value))
}
pval_table <- as.data.frame(pval_table)
colnames(pval_table) <- c("state", "pvalue")
pval_table$pvalue <- as.numeric(as.character(pval_table$pvalue))

saveFile <- paste(outDir, "Supplemental_2_table.txt", sep="")
write.table(pval_table, saveFile, row.names=F, col.names=T, sep="\t", quote=F)





#



#
