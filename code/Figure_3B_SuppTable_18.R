#!/usr/bin/env R
#Figure_3B_SuppTable_18.R

################################################################################
#This script is used to produce Figure 3B as well as Supplemental Table 18.
#
#Run under R 4.1.0
#This script takes as arguments:
# outDir: Path where files are to be written
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# finalAnnotationsTFs: Supplemental Table 1
# abc: Path to ABC connections, provided by Jessie Engreitz
# genomeFile: hg38 genome fasta file.
# meme_dir: a directory containing all meme files generated via provided
#     scripts.  Provided as a zipped file as meme_passed_motifs
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
library(Rsamtools)


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
#This function restricts peaks to the central 101BP around a peak's summit, and
#returns a file formatted for easy conversion into genomicRanges.
################################################################################
formatPeaks <- function(this_peaksFile) {
  this_peaks <- read.table(this_peaksFile, header=F, sep="\t", stringsAsFactors=F)
  this_peaks <- as.data.frame(this_peaks)
  this_peaks[,2] <- as.numeric(as.character(this_peaks[,2]))
  this_peaks[,3] <- as.numeric(as.character(this_peaks[,3]))
  this_peaks[,10] <- as.numeric(as.character(this_peaks[,10]))
  newStarts <- this_peaks[,2]+this_peaks[,10]-50
  newEnds <- this_peaks[,2]+this_peaks[,10]+50

  this_peaks[,2] <- newStarts
  this_peaks[,3] <- newEnds
  colnames(this_peaks)[1:3] <- c("chr", "start", "end")
  return(this_peaks)
}


################################################################################
################################################################################
#Being Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
exprDir <- args[2]
finalAnnotationsTFs <- args[3]
abc <- args[4]
genomeFile <- args[5]
meme_dir <- args[6]

################################################################################
#Identify the different classes-- HOT-HOT, HOT-NOT, NOT-HOT, NOT-NOT
################################################################################


LambertTFs <- read.delim(finalAnnotationsTFs, header=T, sep="\t", stringsAsFactors=F)




abc <- read.delim(abc, header=T, sep="\t", stringsAsFactors=F)
abc <- as.data.frame(abc)
abc <- abc[as.numeric(as.character(abc$distance))>=10000,]

abc[,1] <- as.character(abc[,1])
abc[,2] <- as.numeric(as.character(abc[,2]))
abc[,3] <- as.numeric(as.character(abc[,3]))

abc_gr <- makeGRangesFromDataFrame(abc, ignore.strand=T, keep.extra.columns=F)




TSS_regions <- cbind(abc[,1], abc[,"TargetGeneTSS"], abc[,"TargetGeneTSS"])
colnames(TSS_regions) <- c("chr", "start", "end")
TSS_regions <- as.data.frame(TSS_regions)
TSS_regions[,2] <- as.numeric(as.character(TSS_regions[,2]))
TSS_regions[,3] <- as.numeric(as.character(TSS_regions[,3]))


TSS_regions[,2] <- TSS_regions[,2]-100
TSS_regions[,3] <- TSS_regions[,3]+99
TSS_regions_gr <- makeGRangesFromDataFrame(TSS_regions, ignore.strand=T, keep.extra.columns=T)


peaks <- list.files(exprDir, full.names=T)
peaks <- peaks[grep("Preferred", peaks)]
peaks <- peaks[-grep("H[3|4]K", peaks)]

theTF_names <- c()
enh_peakNum <- rep(0, nrow(abc))
prom_peakNum <- rep(0, nrow(TSS_regions))
enh_TF_names <- c(rep("", nrow(abc)))
prom_TF_names <- c(rep("", nrow(TSS_regions)))
for (i in 1:length(peaks)) {
  theName <- strsplit(peaks[i], split="/")[[1]][length(strsplit(peaks[i], split="/")[[1]])]
  theName <- strsplit(theName, split="_")[[1]][1]
  theName <- gsub("-FLAG", "", theName)
  theName <- gsub("-eGFP", "", theName)
  theTF_names <- c(theTF_names, theName)
  print(paste(i, length(peaks), theName))

  theFile <- formatPeaks(peaks[i])
  theFile_gr <- makeGRangesFromDataFrame(theFile, ignore.strand=T, keep.extra.columns=T)

  promoter_intersect <- as.data.frame(findOverlaps(TSS_regions_gr, theFile_gr))
  prom_peakNum[unique(promoter_intersect[,1])] <- prom_peakNum[unique(promoter_intersect[,1])]+1

  enhancer_intersect <- as.data.frame(findOverlaps(abc_gr, theFile_gr))
  enh_peakNum[unique(enhancer_intersect[,1])] <- enh_peakNum[unique(enhancer_intersect[,1])]+1

  enh_TF_names[unique(enhancer_intersect[,1])] <- paste(enh_TF_names[unique(enhancer_intersect[,1])], theName, sep=",")
  prom_TF_names[unique(promoter_intersect[,1])] <- paste(prom_TF_names[unique(promoter_intersect[,1])], theName, sep=",")


}



annotations <- read.delim(finalAnnotationsTFs, header=T, sep="\t", stringsAsFactors=F)
annotations <- annotations[,c("Factor", "Final.HA.Annotation")]
annotations <- annotations[annotations[,2]=="TF",]


prom_perc_inEnh <- c()
enh_perc_inProm <- c()
jaccard_between <- c()
jaccard_between_TF_only <- c()
for (i in 1:length(enh_TF_names)) {
  enh_Set <- strsplit(enh_TF_names[i], split=",")[[1]]
  enh_Set <- enh_Set[enh_Set!=""]
  prom_Set <- strsplit(prom_TF_names[i], split=",")[[1]]
  prom_Set <- prom_Set[prom_Set!=""]

  prom_perc_inEnh <- c(prom_perc_inEnh, length(prom_Set[prom_Set%in%enh_Set])/length(prom_Set))
  enh_perc_inProm <- c(enh_perc_inProm, length(enh_Set[enh_Set%in%prom_Set])/length(enh_Set))
  jaccard_between <- c(jaccard_between, length(enh_Set[enh_Set%in%prom_Set])/(length(enh_Set[enh_Set%in%prom_Set]) + length(enh_Set[!enh_Set%in%prom_Set]) + length(prom_Set[!prom_Set%in%enh_Set])))

  enh_Set_TFOnly <- enh_Set[enh_Set%in%annotations[,1]]
  prom_Set_TFOnly <- prom_Set[prom_Set%in%annotations[,1]]
  jaccard_between_TF_only <- c(jaccard_between_TF_only, length(enh_Set_TFOnly[enh_Set_TFOnly%in%prom_Set_TFOnly])/(length(enh_Set_TFOnly[enh_Set_TFOnly%in%prom_Set_TFOnly]) + length(enh_Set_TFOnly[!enh_Set_TFOnly%in%prom_Set_TFOnly]) + length(prom_Set_TFOnly[!prom_Set_TFOnly%in%enh_Set_TFOnly])))

}







theCutoff <- floor(length(peaks)/4)


enh_HOT <- rep("NOT_enh", length(enh_peakNum))
enh_HOT[enh_peakNum>theCutoff] <- "HOT_enh"
prom_HOT <- rep("NOT_prom", length(prom_peakNum))
prom_HOT[prom_peakNum>theCutoff] <- "HOT_prom"

data_frame_abc <- cbind(abc[,1:4], abc[,"TargetGene"], abc[,"ABC.Score"], enh_peakNum, prom_peakNum, prom_perc_inEnh, enh_perc_inProm, jaccard_between, jaccard_between_TF_only, enh_HOT, prom_HOT)
colnames(data_frame_abc) <- c("chr", "start", "end", "name", "gene", "ABC.Score", "enh_peakNum", "prom_peakNum", "prom_perc_inEnh", "enh_perc_inProm", "jaccard_between", "jaccard_between_TF_only", "enh_HOT", "prom_HOT")
data_frame_abc <- as.data.frame(data_frame_abc)
data_frame_abc[,1] <- as.character(data_frame_abc[,1])
data_frame_abc[,4] <- as.character(data_frame_abc[,4])
data_frame_abc[,5] <- as.character(data_frame_abc[,5])
data_frame_abc[,13] <- factor(as.character(data_frame_abc[,13]), levels=c("HOT_enh", "NOT_enh"))
data_frame_abc[,14] <- factor(as.character(data_frame_abc[,14]), levels=c("HOT_prom", "NOT_prom"))
for (i in c(2,3,6:12)) { abc[,i] <- as.numeric(as.character(data_frame_abc[,i]))}




HOT <- rep("Neither", nrow(data_frame_abc))
enh_hot <- which(data_frame_abc[,"enh_HOT"]=="HOT_enh")
prom_hot <- which(data_frame_abc[,"prom_HOT"]=="HOT_prom")
both_hot <- enh_hot[enh_hot%in%prom_hot]

HOT[enh_hot] <- "Enhancer"
HOT[prom_hot] <- "Promoter"
HOT[both_hot] <- "Both"

data_frame_abc$HOT <- HOT





################################################################################
#Identify locations of motifs in connections
################################################################################

genome <- FaFile(file=genomeFile)



data_frame_abc_both <- data_frame_abc[data_frame_abc$HOT=="Both",]
TSS_regions_both <- TSS_regions[data_frame_abc$HOT=="Both",]

data_frame_abc_both_gr <- makeGRangesFromDataFrame(data_frame_abc_both[,1:3])
TSS_regions_both_gr <- makeGRangesFromDataFrame(TSS_regions_both[,1:3])

meme_files <- list.files(meme_dir, full.names=T)

motif_number_table_both <- c()

motif_connections_infoTable_both <- c()

for (i in 1:length(peaks)) {
  theName <- strsplit(peaks[i], split="/")[[1]][length(strsplit(peaks[i], split="/")[[1]])]
  this_meme <- gsub(".bed", "", theName)
  this_meme <- gsub(".gz", "", this_meme)
  these_passed_motifs <- meme_files[grep(this_meme, meme_files)]
  if(length(these_passed_motifs)>0 && file.exists(these_passed_motifs)) {
    theMotifs <- read_meme(these_passed_motifs)
    if(length(theMotifs)>0) {

      theName <- strsplit(theName, split="_")[[1]][1]
      theName <- gsub("-FLAG", "", theName)
      theName <- gsub("-eGFP", "", theName)
      print(paste(i, length(peaks), theName))

      theFile <- formatPeaks(peaks[i])
      theFile_gr <- makeGRangesFromDataFrame(theFile, ignore.strand=T, keep.extra.columns=T)

      theSequences <- get_sequence(theFile_gr, genome)
      confirmation_Hits <- runFimo(sequences=theSequences, motifs=theMotifs, meme_path="/cluster/software/meme-5.3.3/bin")

      motifs_intersect <- as.data.frame(findOverlaps(confirmation_Hits, theFile_gr))

      promoter_intersect <- as.data.frame(findOverlaps(TSS_regions_both_gr, theFile_gr))

      enhancer_intersect <- as.data.frame(findOverlaps(data_frame_abc_both_gr, theFile_gr))



      prom_w_peaks <- length(unique(promoter_intersect[,1]))
      enh_w_peaks <- length(unique(enhancer_intersect[,1]))

      prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
      enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

      motif_number_table_both <- rbind(motif_number_table_both, c(theName, length(prom_peaks_w_motif), length(enh_peaks_w_motif)))

      prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
      enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


      both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
      both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
      both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
      both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
      both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

      thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction)

      motif_connections_infoTable_both <- rbind(motif_connections_infoTable_both, thisLine)
    }
  }

}

colnames(motif_connections_infoTable_both) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction")
motif_connections_infoTable_both <- as.data.frame(motif_connections_infoTable_both)
for (i in 3:ncol(motif_connections_infoTable_both)) { motif_connections_infoTable_both[,i] <- as.numeric(as.character(motif_connections_infoTable_both[,i]))}
motif_connections_infoTable_both$logRatio <- log((motif_connections_infoTable_both$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_both$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_both$type <- rep("Both", nrow(motif_connections_infoTable_both))

colnames(motif_number_table_both) <- c("TF", "bothHOT_promoter_motifs_in_peaks", "botHOT_enhancer_motifs_in_peaks")
motif_number_table_both <- as.data.frame(motif_number_table_both)
for (i in 2:3) { motif_number_table_both[,i] <- as.numeric(as.character(motif_number_table_both[,i]))}




################################################################################
#Now for Enhancers
################################################################################


data_frame_abc_enh <- data_frame_abc[data_frame_abc$HOT=="Enhancer",]
TSS_regions_enh <- TSS_regions[data_frame_abc$HOT=="Enhancer",]

data_frame_abc_enh_gr <- makeGRangesFromDataFrame(data_frame_abc_enh[,1:3])
TSS_regions_enh_gr <- makeGRangesFromDataFrame(TSS_regions_enh[,1:3])


motif_number_table_enhancers <- c()

motif_connections_infoTable_enh <- c()

for (i in 1:length(peaks)) {
  theName <- strsplit(peaks[i], split="/")[[1]][length(strsplit(peaks[i], split="/")[[1]])]
  this_meme <- gsub(".bed", "", theName)
  this_meme <- gsub(".gz", "", this_meme)
  these_passed_motifs <- meme_files[grep(this_meme, meme_files)]
  if(length(these_passed_motifs)>0 && file.exists(these_passed_motifs)) {
    theMotifs <- read_meme(these_passed_motifs)
    if(length(theMotifs)>0) {

      theName <- strsplit(theName, split="_")[[1]][1]
      theName <- gsub("-FLAG", "", theName)
      theName <- gsub("-eGFP", "", theName)
      print(paste(i, length(peaks), theName))

      theFile <- formatPeaks(peaks[i])
      theFile_gr <- makeGRangesFromDataFrame(theFile, ignore.strand=T, keep.extra.columns=T)

      theSequences <- get_sequence(theFile_gr, genome)
      confirmation_Hits <- runFimo(sequences=theSequences, motifs=theMotifs, meme_path="/cluster/software/meme-5.3.3/bin")


      motifs_intersect <- as.data.frame(findOverlaps(confirmation_Hits, theFile_gr))

      promoter_intersect <- as.data.frame(findOverlaps(TSS_regions_enh_gr, theFile_gr))

      enhancer_intersect <- as.data.frame(findOverlaps(data_frame_abc_enh_gr, theFile_gr))



      prom_w_peaks <- length(unique(promoter_intersect[,1]))
      enh_w_peaks <- length(unique(enhancer_intersect[,1]))

      prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
      enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

      motif_number_table_enhancers <- rbind(motif_number_table_enhancers, c(theName, length(prom_peaks_w_motif), length(enh_peaks_w_motif)))

      prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
      enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


      both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
      both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
      both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
      both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
      both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

      thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction)

      motif_connections_infoTable_enh <- rbind(motif_connections_infoTable_enh, thisLine)
    }
  }

}



colnames(motif_connections_infoTable_enh) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction")
motif_connections_infoTable_enh <- as.data.frame(motif_connections_infoTable_enh)
for (i in 3:ncol(motif_connections_infoTable_enh)) { motif_connections_infoTable_enh[,i] <- as.numeric(as.character(motif_connections_infoTable_enh[,i]))}
motif_connections_infoTable_enh$logRatio <- log((motif_connections_infoTable_enh$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_enh$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_enh$type <- rep("Enhancer", nrow(motif_connections_infoTable_enh))


colnames(motif_number_table_enhancers) <- c("TF", "enhancerHOT_promoter_motifs_in_peaks", "enhancerHOT_enhancer_motifs_in_peaks")
motif_number_table_enhancers <- as.data.frame(motif_number_table_enhancers)
for (i in 2:3) { motif_number_table_enhancers[,i] <- as.numeric(as.character(motif_number_table_enhancers[,i]))}





################################################################################
#Now for Promoters
################################################################################


data_frame_abc_pro <- data_frame_abc[data_frame_abc$HOT=="Promoter",]
TSS_regions_pro <- TSS_regions[data_frame_abc$HOT=="Promoter",]

data_frame_abc_pro_gr <- makeGRangesFromDataFrame(data_frame_abc_pro[,1:3])
TSS_regions_pro_gr <- makeGRangesFromDataFrame(TSS_regions_pro[,1:3])


motif_number_table_promoters <- c()

motif_connections_infoTable_pro <- c()

for (i in 1:length(peaks)) {
  theName <- strsplit(peaks[i], split="/")[[1]][length(strsplit(peaks[i], split="/")[[1]])]
  this_meme <- gsub(".bed", "", theName)
  this_meme <- gsub(".gz", "", this_meme)
  these_passed_motifs <- meme_files[grep(this_meme, meme_files)]
  if(length(these_passed_motifs)>0 && file.exists(these_passed_motifs)) {
    theMotifs <- read_meme(these_passed_motifs)
    if(length(theMotifs)>0) {

      theName <- strsplit(theName, split="_")[[1]][1]
      theName <- gsub("-FLAG", "", theName)
      theName <- gsub("-eGFP", "", theName)
      print(paste(i, length(peaks), theName))

      theFile <- formatPeaks(peaks[i])
      theFile_gr <- makeGRangesFromDataFrame(theFile, ignore.strand=T, keep.extra.columns=T)

      theSequences <- get_sequence(theFile_gr, genome)
      confirmation_Hits <- runFimo(sequences=theSequences, motifs=theMotifs, meme_path="/cluster/software/meme-5.3.3/bin")


      motifs_intersect <- as.data.frame(findOverlaps(confirmation_Hits, theFile_gr))

      promoter_intersect <- as.data.frame(findOverlaps(TSS_regions_pro_gr, theFile_gr))

      enhancer_intersect <- as.data.frame(findOverlaps(data_frame_abc_pro_gr, theFile_gr))



      prom_w_peaks <- length(unique(promoter_intersect[,1]))
      enh_w_peaks <- length(unique(enhancer_intersect[,1]))

      prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
      enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

      motif_number_table_promoters <- rbind(motif_number_table_promoters, c(theName, length(prom_peaks_w_motif), length(enh_peaks_w_motif)))

      prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
      enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


      both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
      both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
      both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
      both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
      both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

      thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction)

      motif_connections_infoTable_pro <- rbind(motif_connections_infoTable_pro, thisLine)
    }
  }

}



colnames(motif_connections_infoTable_pro) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction")
motif_connections_infoTable_pro <- as.data.frame(motif_connections_infoTable_pro)
for (i in 3:ncol(motif_connections_infoTable_pro)) { motif_connections_infoTable_pro[,i] <- as.numeric(as.character(motif_connections_infoTable_pro[,i]))}
motif_connections_infoTable_pro$logRatio <- log((motif_connections_infoTable_pro$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_pro$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_pro$type <- rep("Promoter", nrow(motif_connections_infoTable_pro))



colnames(motif_number_table_promoters) <- c("TF", "promoterHOT_promoter_motifs_in_peaks", "promoterHOT_enhancer_motifs_in_peaks")
motif_number_table_promoters <- as.data.frame(motif_number_table_promoters)
for (i in 2:3) { motif_number_table_promoters[,i] <- as.numeric(as.character(motif_number_table_promoters[,i]))}



################################################################################
#Now for Neither
################################################################################


data_frame_abc_nei <- data_frame_abc[data_frame_abc$HOT=="Neither",]
TSS_regions_nei <- TSS_regions[data_frame_abc$HOT=="Neither",]

data_frame_abc_nei_gr <- makeGRangesFromDataFrame(data_frame_abc_nei[,1:3])
TSS_regions_nei_gr <- makeGRangesFromDataFrame(TSS_regions_nei[,1:3])


motif_number_table_none <- c()

motif_connections_infoTable_nei <- c()

for (i in 1:length(peaks)) {
  theName <- strsplit(peaks[i], split="/")[[1]][length(strsplit(peaks[i], split="/")[[1]])]
  this_meme <- gsub(".bed", "", theName)
  this_meme <- gsub(".gz", "", this_meme)
  these_passed_motifs <- meme_files[grep(this_meme, meme_files)]
  if(length(these_passed_motifs)>0 && file.exists(these_passed_motifs)) {
    theMotifs <- read_meme(these_passed_motifs)
    if(length(theMotifs)>0) {

      theName <- strsplit(theName, split="_")[[1]][1]
      theName <- gsub("-FLAG", "", theName)
      theName <- gsub("-eGFP", "", theName)
      print(paste(i, length(peaks), theName))

      theFile <- formatPeaks(peaks[i])
      theFile_gr <- makeGRangesFromDataFrame(theFile, ignore.strand=T, keep.extra.columns=T)

      theSequences <- get_sequence(theFile_gr, genome)
      confirmation_Hits <- runFimo(sequences=theSequences, motifs=theMotifs, meme_path="/cluster/software/meme-5.3.3/bin")


      motifs_intersect <- as.data.frame(findOverlaps(confirmation_Hits, theFile_gr))

      promoter_intersect <- as.data.frame(findOverlaps(TSS_regions_nei_gr, theFile_gr))

      enhancer_intersect <- as.data.frame(findOverlaps(data_frame_abc_nei_gr, theFile_gr))



      prom_w_peaks <- length(unique(promoter_intersect[,1]))
      enh_w_peaks <- length(unique(enhancer_intersect[,1]))

      prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
      enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

      motif_number_table_none <- rbind(motif_number_table_none, c(theName, length(prom_peaks_w_motif), length(enh_peaks_w_motif)))

      prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
      enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


      both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
      both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
      both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
      both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
      both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

      thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction)

      motif_connections_infoTable_nei <- rbind(motif_connections_infoTable_nei, thisLine)
    }
  }

}



colnames(motif_connections_infoTable_nei) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction")
motif_connections_infoTable_nei <- as.data.frame(motif_connections_infoTable_nei)
for (i in 3:ncol(motif_connections_infoTable_nei)) { motif_connections_infoTable_nei[,i] <- as.numeric(as.character(motif_connections_infoTable_nei[,i]))}
motif_connections_infoTable_nei$logRatio <- log((motif_connections_infoTable_nei$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_nei$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_nei$type <- rep("Neither", nrow(motif_connections_infoTable_nei))



colnames(motif_number_table_none) <- c("TF", "neitherHOT_promoter_motifs_in_peaks", "neitherHOT_enhancer_motifs_in_peaks")
motif_number_table_none <- as.data.frame(motif_number_table_none)
for (i in 2:3) { motif_number_table_none[,i] <- as.numeric(as.character(motif_number_table_none[,i]))}




################################################################################
#Without regard to where HOTs are...
################################################################################



data_frame_abc_gr <- makeGRangesFromDataFrame(data_frame_abc[,1:3])
TSS_regions_gr <- makeGRangesFromDataFrame(TSS_regions[,1:3])


motif_number_table_all <- c()

motif_connections_infoTable_all <- c()

for (i in 1:length(peaks)) {
  theName <- strsplit(peaks[i], split="/")[[1]][length(strsplit(peaks[i], split="/")[[1]])]
  this_meme <- gsub(".bed", "", theName)
  this_meme <- gsub(".gz", "", this_meme)
  these_passed_motifs <- meme_files[grep(this_meme, meme_files)]
  if(length(these_passed_motifs)>0 && file.exists(these_passed_motifs)) {
    theMotifs <- read_meme(these_passed_motifs)
    if(length(theMotifs)>0) {

      theName <- strsplit(theName, split="_")[[1]][1]
      theName <- gsub("-FLAG", "", theName)
      theName <- gsub("-eGFP", "", theName)
      print(paste(i, length(peaks), theName))

      theFile <- formatPeaks(peaks[i])
      theFile_gr <- makeGRangesFromDataFrame(theFile, ignore.strand=T, keep.extra.columns=T)

      theSequences <- get_sequence(theFile_gr, genome)
      confirmation_Hits <- runFimo(sequences=theSequences, motifs=theMotifs, meme_path="/cluster/software/meme-5.3.3/bin")


      motifs_intersect <- as.data.frame(findOverlaps(confirmation_Hits, theFile_gr))

      promoter_intersect <- as.data.frame(findOverlaps(TSS_regions_gr, theFile_gr))

      enhancer_intersect <- as.data.frame(findOverlaps(data_frame_abc_gr, theFile_gr))



      prom_w_peaks <- length(unique(promoter_intersect[,1]))
      enh_w_peaks <- length(unique(enhancer_intersect[,1]))

      prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
      enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

      motif_number_table_all <- rbind(motif_number_table_all, c(theName, length(prom_peaks_w_motif), length(enh_peaks_w_motif)))

      prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
      enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


      both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
      both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
      both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
      both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
      both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

      thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction)

      motif_connections_infoTable_all <- rbind(motif_connections_infoTable_all, thisLine)
    }
  }

}



colnames(motif_connections_infoTable_all) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction")
motif_connections_infoTable_all <- as.data.frame(motif_connections_infoTable_all)
for (i in 3:ncol(motif_connections_infoTable_all)) { motif_connections_infoTable_all[,i] <- as.numeric(as.character(motif_connections_infoTable_all[,i]))}
motif_connections_infoTable_all$logRatio <- log((motif_connections_infoTable_all$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_all$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_all$type <- rep("NoRegard", nrow(motif_connections_infoTable_all))



colnames(motif_number_table_all) <- c("TF", "irrespectiveHOT_promoter_motifs_in_peaks", "irrespectiveHOT_enhancer_motifs_in_peaks")
motif_number_table_all <- as.data.frame(motif_number_table_all)
for (i in 2:3) { motif_number_table_all[,i] <- as.numeric(as.character(motif_number_table_all[,i]))}





################################################################################
#Make a shared data frame for all of them.
################################################################################

motifs_in_peaks_count_table <- cbind(motif_number_table_both, motif_number_table_promoters[,2:3], motif_number_table_enhancers[,2:3], motif_number_table_none[,2:3], motif_number_table_all[,2:3])

motifs_in_peaks_count_table_TFsOnly <- motifs_in_peaks_count_table[motifs_in_peaks_count_table[,1]%in%LambertTFs[LambertTFs[,5]=="TF",2],]

saveFile <- paste(outDir, "Supplemental_Table_18.txt", sep="")

write.table(motifs_in_peaks_count_table_TFsOnly, saveFile, row.names=F, col.names=T, sep="\t", quote=F)


################################################################################
#Make a shared data frame for all of them.
################################################################################

motif_connections_infoTable_all <- rbind(motif_connections_infoTable_both, motif_connections_infoTable_enh, motif_connections_infoTable_pro, motif_connections_infoTable_nei)

colnames(motif_connections_infoTable_all)[11] <- "HOT"


################################################################################
#Make a figure of motif distributions, restricted to cases of at least
#100 peaks to avoid cases of small numbers making things appear skewed
#in an outsized way.
################################################################################

motif_connections_infoTable_all_many <- motif_connections_infoTable_all[motif_connections_infoTable_all$cons_w_dual_peaks>=100,]

typeVec <- unique(motif_connections_infoTable_all$HOT)
custom_col <- c("blue", "orange", "grey", "red")

saveFile <- paste(outDir, "Moyers_Figure3B.pdf", sep="")
p <- ggplot(motif_connections_infoTable_all_many, aes(x=logRatio, color=HOT)) + geom_density(size=2) + blank_theme +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=20), legend.text=element_text(size=15)) + xlab("Motif in Enhancer        Motif in Promoter") + xlim(-2,2) +
  scale_color_manual(name="Category", values = custom_col, na.value="grey50")
ggsave(saveFile)

pval_table <- matrix(nrow=4, ncol=4)
colnames(pval_table) <- typeVec
rownames(pval_table) <- typeVec
for (i in 1:(length(typeVec)-1)) {
  for (j in (i+1):length(typeVec)) {
    vec1 <- motif_connections_infoTable_all_many[motif_connections_infoTable_all_many$HOT==typeVec[i],"logRatio"]
    vec2 <- motif_connections_infoTable_all_many[motif_connections_infoTable_all_many$HOT==typeVec[j],"logRatio"]
    theTest <- ks.test(vec1, vec2)
    pval_table[typeVec[i], typeVec[j]] <- theTest$p.value
  }
}
saveFile <- paste(outDir, "Supplemental_Table_19_3B.txt", sep="")
write.table(pval_table, saveFile, row.names=T, col.names=T, sep="\t", quote=F)


motif_connections_infoTable_all_many$Ratio <- exp(motif_connections_infoTable_all_many$logRatio)

saveFile <- paste(outDir, "Moyers_Figure3B_ratio.pdf", sep="")
p <- ggplot(motif_connections_infoTable_all_many, aes(x=Ratio, color=HOT)) + geom_density(size=2) + blank_theme +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=20), legend.text=element_text(size=15)) + xlab("Motif in Enhancer        Motif in Promoter") + xlim(0,3) +
  scale_color_manual(name="Category", values = custom_col, na.value="grey50")
ggsave(saveFile)


#
