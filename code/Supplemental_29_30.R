#!/usr/bin/env R
#Supplemental_29_30.R

################################################################################
#This script is used to produce Supplemental Figures 29 and 30
#Note that, unlike most other figure-producing scripts, this one does
#not specifically name the figure "Supplemental_XX.pdf", but produces a filename
#that is descriptive of its contents.
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
# cgi_promoters: an annotation of promoters' cgi status, provided as a file
#     as proms_v103_annotated_GR_like.tsv
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
cgi_promoters <- args[7]

################################################################################
#Identify the different classes-- HOT-HOT, HOT-NOT, NOT-HOT, NOT-NOT
################################################################################

LambertTFs <- read.table(finalAnnotationsTFs, header=T, sep="\t", stringsAsFactors=F)


abc <- read.delim(abc, header=T, sep="\t", stringsAsFactors=F)
abc <- as.data.frame(abc)
abc <- abc[as.numeric(as.character(abc$distance))>=10000,]

abc[,1] <- as.character(abc[,1])
abc[,2] <- as.numeric(as.character(abc[,2]))
abc[,3] <- as.numeric(as.character(abc[,3]))

abc_gr <- makeGRangesFromDataFrame(abc, ignore.strand=T, keep.extra.columns=F)



cgi_promoters <- read.table(cgi_promoters, header=T, sep="\t", stringsAsFactors=F)

cgi_promoters_CGI_bivalent <- cgi_promoters[cgi_promoters$CGI_class=="CGI_bivalent",]
cgi_promoters_CGI_bivalent <- cgi_promoters_CGI_bivalent[,c("seqnames", "start", "end", "CGI_class")]
cgi_promoters_CGI_constit <- cgi_promoters[cgi_promoters$CGI_class=="CGI_constit",]
cgi_promoters_CGI_constit <- cgi_promoters_CGI_constit[,c("seqnames", "start", "end", "CGI_class")]
colnames(cgi_promoters_CGI_bivalent) <- c("chr", "start", "end", "type")
colnames(cgi_promoters_CGI_constit) <- c("chr", "start", "end", "type")
cgi_promoters_CGI_bivalent_gr <- makeGRangesFromDataFrame(cgi_promoters_CGI_bivalent)
cgi_promoters_CGI_constit_gr <- makeGRangesFromDataFrame(cgi_promoters_CGI_constit)



theType <- rep("nonCGI", nrow(abc))
bivalentIntersect <- as.data.frame(findOverlaps(abc_gr, cgi_promoters_CGI_bivalent_gr))
constitutiveIntersect <- as.data.frame(findOverlaps(abc_gr, cgi_promoters_CGI_constit_gr))
theType[unique(bivalentIntersect[,1])] <- "bivalentCGI"
theType[unique(constitutiveIntersect[,1])] <- "constitCGI"
abc$cgiType <- theType


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

data_frame_abc <- cbind(abc[,1:4], abc[,"TargetGene"], abc[,"ABC.Score"], abc[,"cgiType"], enh_peakNum, prom_peakNum, prom_perc_inEnh, enh_perc_inProm, jaccard_between, jaccard_between_TF_only, enh_HOT, prom_HOT)
colnames(data_frame_abc) <- c("chr", "start", "end", "name", "gene", "ABC.Score", "cgiType", "enh_peakNum", "prom_peakNum", "prom_perc_inEnh", "enh_perc_inProm", "jaccard_between", "jaccard_between_TF_only", "enh_HOT", "prom_HOT")
data_frame_abc <- as.data.frame(data_frame_abc)
data_frame_abc[,1] <- as.character(data_frame_abc[,1])
data_frame_abc[,4] <- as.character(data_frame_abc[,4])
data_frame_abc[,5] <- as.character(data_frame_abc[,5])
data_frame_abc[,6] <- as.character(data_frame_abc[,6])
data_frame_abc[,14] <- factor(as.character(data_frame_abc[,14]), levels=c("HOT_enh", "NOT_enh"))
data_frame_abc[,15] <- factor(as.character(data_frame_abc[,15]), levels=c("HOT_prom", "NOT_prom"))
for (i in c(2,3,7:13)) { abc[,i] <- as.numeric(as.character(data_frame_abc[,i]))}



HOT <- rep("Neither", nrow(data_frame_abc))
enh_hot <- which(data_frame_abc[,"enh_HOT"]=="HOT_enh")
prom_hot <- which(data_frame_abc[,"prom_HOT"]=="HOT_prom")
both_hot <- enh_hot[enh_hot%in%prom_hot]

HOT[enh_hot] <- "Enhancer"
HOT[prom_hot] <- "Promoter"
HOT[both_hot] <- "Both"

data_frame_abc$HOT <- HOT



TSS_regions$HOT <- data_frame_abc$HOT
TSS_regions$cgiType <- data_frame_abc$cgiType


################################################################################
#Identify locations of motifs in connections
################################################################################

cgiTypes <- c("nonCGI", "bivalentCGI", "constitutiveCGI")


genome <- FaFile(file=genomeFile)



data_frame_abc_both <- data_frame_abc[data_frame_abc$HOT=="Both",]

data_frame_abc_both_nonCGI <- data_frame_abc_both[data_frame_abc_both$cgiType=="nonCGI",]
data_frame_abc_both_bivalentCGI <- data_frame_abc_both[data_frame_abc_both$cgiType=="bivalentCGI",]
data_frame_abc_both_constitCGI <- data_frame_abc_both[data_frame_abc_both$cgiType=="constitCGI",]


TSS_regions_both <- TSS_regions[data_frame_abc$HOT=="Both",]

TSS_regions_both_nonCGI <- TSS_regions_both[TSS_regions_both$cgiType=="nonCGI",]
TSS_regions_both_bivalentCGI <- TSS_regions_both[TSS_regions_both$cgiType=="bivalentCGI",]
TSS_regions_both_constitCGI <- TSS_regions_both[TSS_regions_both$cgiType=="constitCGI",]



data_frame_abc_both_nonCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_both_nonCGI[,1:3])
data_frame_abc_both_bivalentCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_both_bivalentCGI[,1:3])
data_frame_abc_both_constitCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_both_constitCGI[,1:3])

TSS_regions_both_gr <- makeGRangesFromDataFrame(TSS_regions_both[,1:3])
TSS_regions_both_nonCGI_gr <- makeGRangesFromDataFrame(TSS_regions_both_nonCGI[,1:3])
TSS_regions_both_bivalentCGI_gr <- makeGRangesFromDataFrame(TSS_regions_both_bivalentCGI[,1:3])
TSS_regions_both_constitCGI_gr <- makeGRangesFromDataFrame(TSS_regions_both_constitCGI[,1:3])


meme_files <- list.files(meme_dir, full.names=T)

promList <- list(TSS_regions_both_nonCGI_gr, TSS_regions_both_bivalentCGI_gr, TSS_regions_both_constitCGI_gr)
enhList <- list(data_frame_abc_both_nonCGI_gr, data_frame_abc_both_bivalentCGI_gr, data_frame_abc_both_constitCGI_gr)

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

      for (j in 1:length(cgiTypes)) {

        promoter_intersect <- as.data.frame(findOverlaps(promList[[j]], theFile_gr))

        enhancer_intersect <- as.data.frame(findOverlaps(enhList[[j]], theFile_gr))

        prom_w_peaks <- length(unique(promoter_intersect[,1]))
        enh_w_peaks <- length(unique(enhancer_intersect[,1]))

        prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
        enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

        prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
        enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


        both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
        both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
        both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
        both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
        both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

        thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction, cgiTypes[j])

        motif_connections_infoTable_both <- rbind(motif_connections_infoTable_both, thisLine)
      }
    }
  }

}

colnames(motif_connections_infoTable_both) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction", "cgiType")
motif_connections_infoTable_both <- as.data.frame(motif_connections_infoTable_both)
for (i in 3:(ncol(motif_connections_infoTable_both)-1)) { motif_connections_infoTable_both[,i] <- as.numeric(as.character(motif_connections_infoTable_both[,i]))}
motif_connections_infoTable_both$logRatio <- log((motif_connections_infoTable_both$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_both$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_both$type <- rep("Both", nrow(motif_connections_infoTable_both))





################################################################################
#Now for Enhancers
################################################################################


data_frame_abc_enh <- data_frame_abc[data_frame_abc$HOT=="Enhancer",]

data_frame_abc_enh_nonCGI <- data_frame_abc_enh[data_frame_abc_enh$cgiType=="nonCGI",]
data_frame_abc_enh_bivalentCGI <- data_frame_abc_enh[data_frame_abc_enh$cgiType=="bivalentCGI",]
data_frame_abc_enh_constitCGI <- data_frame_abc_enh[data_frame_abc_enh$cgiType=="constitCGI",]

TSS_regions_enh <- TSS_regions[data_frame_abc$HOT=="Enhancer",]

TSS_regions_enh_nonCGI <- TSS_regions_enh[TSS_regions_enh$cgiType=="nonCGI",]
TSS_regions_enh_bivalentCGI <- TSS_regions_enh[TSS_regions_enh$cgiType=="bivalentCGI",]
TSS_regions_enh_constitCGI <- TSS_regions_enh[TSS_regions_enh$cgiType=="constitCGI",]

data_frame_abc_enh_gr <- makeGRangesFromDataFrame(data_frame_abc_enh[,1:3])

data_frame_abc_enh_nonCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_enh_nonCGI[,1:3])
data_frame_abc_enh_bivalentCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_enh_bivalentCGI[,1:3])
data_frame_abc_enh_constitCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_enh_constitCGI[,1:3])


TSS_regions_enh_nonCGI_gr <- makeGRangesFromDataFrame(TSS_regions_enh_nonCGI[,1:3])
TSS_regions_enh_bivalentCGI_gr <- makeGRangesFromDataFrame(TSS_regions_enh_bivalentCGI[,1:3])
TSS_regions_enh_constitCGI_gr <- makeGRangesFromDataFrame(TSS_regions_enh_constitCGI[,1:3])



promList <- list(TSS_regions_enh_nonCGI_gr, TSS_regions_enh_bivalentCGI_gr, TSS_regions_enh_constitCGI_gr)
enhList <- list(data_frame_abc_enh_nonCGI_gr, data_frame_abc_enh_bivalentCGI_gr, data_frame_abc_enh_constitCGI_gr)


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

      for (j in 1:length(cgiTypes)) {

        promoter_intersect <- as.data.frame(findOverlaps(promList[[j]], theFile_gr))

        enhancer_intersect <- as.data.frame(findOverlaps(enhList[[j]], theFile_gr))

        prom_w_peaks <- length(unique(promoter_intersect[,1]))
        enh_w_peaks <- length(unique(enhancer_intersect[,1]))

        prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
        enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

        prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
        enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


        both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
        both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
        both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
        both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
        both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

        thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction, cgiTypes[j])

        motif_connections_infoTable_enh <- rbind(motif_connections_infoTable_enh, thisLine)
      }
    }
  }

}



colnames(motif_connections_infoTable_enh) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction", "cgiType")
motif_connections_infoTable_enh <- as.data.frame(motif_connections_infoTable_enh)
for (i in 3:(ncol(motif_connections_infoTable_enh)-1)) { motif_connections_infoTable_enh[,i] <- as.numeric(as.character(motif_connections_infoTable_enh[,i]))}
motif_connections_infoTable_enh$logRatio <- log((motif_connections_infoTable_enh$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_enh$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_enh$type <- rep("Enhancer", nrow(motif_connections_infoTable_enh))



################################################################################
#Now for Promoters
################################################################################


data_frame_abc_pro <- data_frame_abc[data_frame_abc$HOT=="Promoter",]
data_frame_abc_pro_nonCGI <- data_frame_abc_pro[data_frame_abc_pro$cgiType=="nonCGI",]
data_frame_abc_pro_bivalentCGI <- data_frame_abc_pro[data_frame_abc_pro$cgiType=="bivalentCGI",]
data_frame_abc_pro_constitCGI <- data_frame_abc_pro[data_frame_abc_pro$cgiType=="constitCGI",]

TSS_regions_pro <- TSS_regions[data_frame_abc$HOT=="Promoter",]
TSS_regions_pro_nonCGI <- TSS_regions_pro[TSS_regions_pro$cgiType=="nonCGI",]
TSS_regions_pro_bivalentCGI <- TSS_regions_pro[TSS_regions_pro$cgiType=="bivalentCGI",]
TSS_regions_pro_constitCGI <- TSS_regions_pro[TSS_regions_pro$cgiType=="constitCGI",]


data_frame_abc_pro_nonCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_pro_nonCGI[,1:3])
data_frame_abc_pro_bivalentCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_pro_bivalentCGI[,1:3])
data_frame_abc_pro_constitCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_pro_constitCGI[,1:3])


TSS_regions_pro_nonCGI_gr <- makeGRangesFromDataFrame(TSS_regions_pro_nonCGI[,1:3])
TSS_regions_pro_bivalentCGI_gr <- makeGRangesFromDataFrame(TSS_regions_pro_bivalentCGI[,1:3])
TSS_regions_pro_constitCGI_gr <- makeGRangesFromDataFrame(TSS_regions_pro_constitCGI[,1:3])



promList <- list(TSS_regions_pro_nonCGI_gr, TSS_regions_pro_bivalentCGI_gr, TSS_regions_pro_constitCGI_gr)
enhList <- list(data_frame_abc_pro_nonCGI_gr, data_frame_abc_pro_bivalentCGI_gr, data_frame_abc_pro_constitCGI_gr)


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

      for (j in 1:length(cgiTypes)) {

        promoter_intersect <- as.data.frame(findOverlaps(promList[[j]], theFile_gr))

        enhancer_intersect <- as.data.frame(findOverlaps(enhList[[j]], theFile_gr))

        prom_w_peaks <- length(unique(promoter_intersect[,1]))
        enh_w_peaks <- length(unique(enhancer_intersect[,1]))

        prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
        enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

        prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
        enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


        both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
        both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
        both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
        both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
        both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

        thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction, cgiTypes[j])

        motif_connections_infoTable_pro <- rbind(motif_connections_infoTable_pro, thisLine)
      }
    }
  }

}



colnames(motif_connections_infoTable_pro) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction", "cgiType")
motif_connections_infoTable_pro <- as.data.frame(motif_connections_infoTable_pro)
for (i in 3:(ncol(motif_connections_infoTable_pro)-1)) { motif_connections_infoTable_pro[,i] <- as.numeric(as.character(motif_connections_infoTable_pro[,i]))}
motif_connections_infoTable_pro$logRatio <- log((motif_connections_infoTable_pro$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_pro$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_pro$type <- rep("Promoter", nrow(motif_connections_infoTable_pro))





################################################################################
#Now for Neither
################################################################################


data_frame_abc_nei <- data_frame_abc[data_frame_abc$HOT=="Neither",]
data_frame_abc_nei_nonCGI <- data_frame_abc_nei[data_frame_abc_nei$cgiType=="nonCGI",]
data_frame_abc_nei_bivalentCGI <- data_frame_abc_nei[data_frame_abc_nei$cgiType=="bivalentCGI",]
data_frame_abc_nei_constitCGI <- data_frame_abc_nei[data_frame_abc_nei$cgiType=="constitCGI",]

TSS_regions_nei <- TSS_regions[data_frame_abc$HOT=="Neither",]
TSS_regions_nei_nonCGI <- TSS_regions_nei[TSS_regions_nei$cgiType=="nonCGI",]
TSS_regions_nei_bivalentCGI <- TSS_regions_nei[TSS_regions_nei$cgiType=="bivalentCGI",]
TSS_regions_nei_constitCGI <- TSS_regions_nei[TSS_regions_nei$cgiType=="constitCGI",]


data_frame_abc_nei_nonCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_nei_nonCGI[,1:3])
data_frame_abc_nei_bivalentCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_nei_bivalentCGI[,1:3])
data_frame_abc_nei_constitCGI_gr <- makeGRangesFromDataFrame(data_frame_abc_nei_constitCGI[,1:3])


TSS_regions_nei_nonCGI_gr <- makeGRangesFromDataFrame(TSS_regions_nei_nonCGI[,1:3])
TSS_regions_nei_bivalentCGI_gr <- makeGRangesFromDataFrame(TSS_regions_nei_bivalentCGI[,1:3])
TSS_regions_nei_constitCGI_gr <- makeGRangesFromDataFrame(TSS_regions_nei_constitCGI[,1:3])


promList <- list(TSS_regions_nei_nonCGI_gr, TSS_regions_nei_bivalentCGI_gr, TSS_regions_nei_constitCGI_gr)
enhList <- list(data_frame_abc_nei_nonCGI_gr, data_frame_abc_nei_bivalentCGI_gr, data_frame_abc_nei_constitCGI_gr)



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

      for (j in 1:length(cgiTypes)) {

        promoter_intersect <- as.data.frame(findOverlaps(promList[[j]], theFile_gr))

        enhancer_intersect <- as.data.frame(findOverlaps(enhList[[j]], theFile_gr))

        prom_w_peaks <- length(unique(promoter_intersect[,1]))
        enh_w_peaks <- length(unique(enhancer_intersect[,1]))

        prom_peaks_w_motif <- promoter_intersect[promoter_intersect[,2]%in%motifs_intersect[,2],2]
        enh_peaks_w_motif <- enhancer_intersect[enhancer_intersect[,2]%in%motifs_intersect[,2],2]

        prom_peaks_w_motif_fraction <- length(unique(prom_peaks_w_motif))/prom_w_peaks
        enh_peaks_w_motif_fraction <- length(unique(enh_peaks_w_motif))/enh_w_peaks


        both_w_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],]
        both_enh_peaks <- enhancer_intersect[enhancer_intersect[,1]%in%promoter_intersect[,1],2]
        both_prom_peaks <- promoter_intersect[promoter_intersect[,1]%in%enhancer_intersect[,1],2]
        both_enh_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_enh_peaks,2])) / length(unique(both_enh_peaks))
        both_prom_peaks_fraction <- length(unique(motifs_intersect[motifs_intersect[,2]%in%both_prom_peaks,2])) / length(unique(both_prom_peaks))

        thisLine <- c(theName, peaks[i], prom_w_peaks, enh_w_peaks, prom_peaks_w_motif_fraction, enh_peaks_w_motif_fraction, length(unique(both_w_peaks[,1])), both_prom_peaks_fraction, both_enh_peaks_fraction, cgiTypes[j])

        motif_connections_infoTable_nei <- rbind(motif_connections_infoTable_nei, thisLine)
      }
    }
  }

}



colnames(motif_connections_infoTable_nei) <- c("TF", "expr", "proms_w_peaks", "enh_w_peaks", "prom_peaks_w_motif_fraction", "enh_peaks_w_motif_fraction", "cons_w_dual_peaks", "both_prom_peaks_fraction", "both_enh_peaks_fraction", "cgiType")
motif_connections_infoTable_nei <- as.data.frame(motif_connections_infoTable_nei)
for (i in 3:(ncol(motif_connections_infoTable_nei)-1)) { motif_connections_infoTable_nei[,i] <- as.numeric(as.character(motif_connections_infoTable_nei[,i]))}
motif_connections_infoTable_nei$logRatio <- log((motif_connections_infoTable_nei$both_prom_peaks_fraction + 0.1) / (motif_connections_infoTable_nei$both_enh_peaks_fraction + 0.1))

motif_connections_infoTable_nei$type <- rep("Neither", nrow(motif_connections_infoTable_nei))



################################################################################
#Make a shared data frame for all of them.
################################################################################

motif_connections_infoTable_all <- rbind(motif_connections_infoTable_both, motif_connections_infoTable_enh, motif_connections_infoTable_pro, motif_connections_infoTable_nei)

colnames(motif_connections_infoTable_all)[12] <- "HOT"



################################################################################
#Make a figure of motif distributions, restricted to cases of at least
#100 peaks to avoid cases of small numbers making things appear skewed
#in an outsized way.
################################################################################



typeVec <- c("Both", "Promoter", "Enhancer", "Neither")

custom_col <- c("blue", "yellow", "grey", "red")


for (k in 1:length(cgiTypes)) {
  this_motif_connections_infoTable_all <- motif_connections_infoTable_all[motif_connections_infoTable_all$cgiType==cgiTypes[k],]
  if(nrow(this_motif_connections_infoTable_all)>=3) {
    saveFile <- paste(outDir, "Distal_All_logRatio_promVenh_cgiAware_", cgiTypes[k], "_TFsOnly.pdf", sep="")
    p <- ggplot(this_motif_connections_infoTable_all, aes(x=logRatio, color=HOT)) + geom_density(size=2) + blank_theme +
      theme(axis.text=element_text(size=20), axis.title=element_text(size=25), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=20), legend.text=element_text(size=15)) + xlab("Motif in Enhancer        Motif in Promoter") + xlim(-2,2) +
      scale_color_manual(name="Category", values = custom_col, na.value="grey50")
    ggsave(saveFile)

    pval_table <- matrix(nrow=4, ncol=4)
    colnames(pval_table) <- typeVec
    rownames(pval_table) <- typeVec
    for (i in 1:(length(typeVec)-1)) {
      for (j in (i+1):length(typeVec)) {
        vec1 <- this_motif_connections_infoTable_all[this_motif_connections_infoTable_all$HOT==typeVec[i],"logRatio"]
        vec2 <- this_motif_connections_infoTable_all[this_motif_connections_infoTable_all$HOT==typeVec[j],"logRatio"]
        if(length(vec1)>=2 && length(vec2)>=3) {
          theTest <- ks.test(vec1, vec2)
          pval_table[typeVec[i], typeVec[j]] <- theTest$p.value
        }
      }
    }
    saveFile <- paste(outDir, "Supplemental_Table_19", cgiTypes[k], ".txt", sep="")
    write.table(pval_table, saveFile, row.names=T, col.names=T, sep="\t", quote=F)
    print((cgiTypes[k]))
    print(pval_table)
  }

}





#
