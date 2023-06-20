#Fig4_C_D_E.R

################################################################################
#Run under R/4.1.1
#This script takes as arguments:
# outDir: Path where files are to be written
# finalAnnotationsTFs: Supplemental Table 1
# cCREs: Path to the Version 4 cCREs provided by Jill Moore in June 2022.
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# thePromoters: file containing refseq TSS with gene names +/- 1kb,
#     provided for download as refseq_genes_unique_TSS.bed
# expr1 and expr2: Path to expression levels of genes in HepG2,
#     accession numbers ENCFF533XPJ and ENCFF321JIT
# hotSites: path to a bed-like file of HOT sites, provided as HOT_Sites.bed
# genomeFile: hg38 genome fasta file.
# hiccup: path to HiC data found under accession ENCFF264RQT
#
################################################################################


################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(ggplot2)
library(rtracklayer)
library(universalmotif)
library(memes)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(plyranges)
library("biomaRt")

options("scipen"=100)




################################################################################
################################################################################
#Define Functions
################################################################################
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


running_meme_pipeline <- function(thePeaks, hotSites_gr, genome, this_outDir) {
  #thePeaks <- sMAF_Cofactor_peaks

  #Ensure the outdir exists.
  if(!file.exists(this_outDir)) { system(paste("mkdir ", this_outDir, sep=""))}

  #Remove HOT sites.
  thePeaks_gr <- makeGRangesFromDataFrame(thePeaks)
  theOverlap <- as.data.frame(findOverlaps(thePeaks_gr, hotSites_gr))
  thePeaks <- thePeaks[-unique(theOverlap[,1]),]

  #Sort by Peak Score. Take the first 500 peaks.
  thePeaks <- thePeaks[order(thePeaks[,7], decreasing=T),]
  thePeaks_top <- thePeaks[1:500,]
  thePeaks_top <- thePeaks_top[!is.na(thePeaks_top[,1]),]
  thePeaks_top_gr <- makeGRangesFromDataFrame(thePeaks_top)

  #Get control sequences.
  saveFile <- paste(this_outDir, "top_500_peaks.bed", sep="")
  write.table(thePeaks_top, saveFile, row.names=F, col.names=F, sep="\t", quote=F)
  control_file <- paste(this_outDir, "top_500_peaks_control.bed", sep="")
  theCommand <- paste("module load g/python/2.7.15; module load g/python/2.7-modules; python /gpfs/gpfs1/home/bmoyers/Biotrain_2019/Scripts/nullseq_generate.py -x 1 -m 1000 -r 1 -o ", control_file, " ", saveFile, " hg38 /cluster/lab/myers/reference/nullseq_indices/hg38_indices/", sep="")
  system(theCommand)

  controlPeaks <- read.table(control_file, header=F, sep="\t", stringsAsFactors=F)
  controlPeaks <- as.data.frame(controlPeaks)
  colnames(controlPeaks)[1:3] <- c("chr", "start", "end")
  controlPeaks_gr <- makeGRangesFromDataFrame(controlPeaks)

  #get a BioStrings object for each of the experimental and control sets.
  inputSeqs <- get_sequence(thePeaks_top_gr, genome)
  controlSeqs <- get_sequence(controlPeaks_gr, genome)

  meme_outDir <- paste(this_outDir, "meme", sep="")
  system("module load g/meme/5.3.3")
  theMotifs <- runMeme(input=inputSeqs, control = controlSeqs, outdir = meme_outDir, alph = "dna",
    silent = FALSE, meme_path = "/cluster/software/meme-5.3.3/bin", nmotifs=5, minw=6, maxw=50, objfun="de")


}


################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
finalAnnotationsTFs <- args[2]
cCREs <- args[3]
exprDir <- args[4]
thePromoters <- args[5]
expr1 <- args[6]
expr2 <- args[7]
hotSites <- args[8]
genomeFile <- args[9]
hiccup <- args[10]


#Identify all of the peaks.
#pull out the sMAFs.
#pull out the NFE2 experiment.
#restrict to TFs only, which will require lambert table and replacing CSDA/YBX3.



LambertTFs <- read.table(finalAnnotationsTFs, header=T, sep="\t", stringsAsFactors=F)
LambertTFs <- LambertTFs[LambertTFs[,3]=="Preferred",]
LambertTFs <- LambertTFs[LambertTFs[,5]=="TF",]


#Read in the cCREs for later use.
cCREs <- read.table(cCREs, header=F, sep="\t", stringsAsFactors=F)
cCREs <- cCREs[,c(1:3,6)]
colnames(cCREs)[1:4] <- c("chr", "start", "end", "state")
cCREs <- cCREs[cCREs[,4]!="TF",]
cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=TRUE, keep.extra.columns=F)





peakFiles <- list.files(exprDir, full.names=T)
peakFiles <- peakFiles[grep("Preferred", peakFiles)]

TFs_table <- matrix(nrow=length(peakFiles), ncol=2)
for (i in 1:length(peakFiles)) {
  thisTF <- strsplit(peakFiles[i], split="/")[[1]]
  thisTF <- strsplit(thisTF[length(thisTF)], split="_")[[1]][1]
  thisTF <- gsub("-FLAG", "", thisTF)
  thisTF <- gsub("-eGFP", "", thisTF)
  TFs_table[i,1] <- thisTF
  TFs_table[i,2] <- peakFiles[i]
}


dim(TFs_table[TFs_table[,1]%in%LambertTFs[LambertTFs[,4]=="Yes",2],])

TFs_table <- TFs_table[TFs_table[,1]%in%LambertTFs[LambertTFs[,4]=="Yes",2],]

TFs_table <- as.data.frame(TFs_table)
colnames(TFs_table) <- c("TF", "PeakPath")

sMAF_table <- TFs_table[TFs_table[,1]%in%c("MAFF", "MAFK", "MAFG"),]
################################################################################
#Cofactors ID'd in fig 5 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911266/
################################################################################
Cofactor_table <- TFs_table[TFs_table[,1]%in%c("NFE2", "NRF1", "NRF2", "NRF3", "BACH1", "BACH2", "JUN", "FOS"),]
TFs_table <- TFs_table[!TFs_table[,1]%in%c("MAFF", "MAFK", "MAFG", "NFE2", "NRF1", "NRF2", "NRF3", "BACH1", "BACH2", "JUN", "FOS"),]


sMAF_peaks <- c()
for (i in 1:nrow(sMAF_table)) {
  thesePeaks <- read.table(sMAF_table[i,2], header=F, sep="\t", stringsAsFactors=F)
  thesePeaks <- cbind(thesePeaks, rep(sMAF_table[i,1], nrow(thesePeaks)))
  colnames(thesePeaks)[ncol(thesePeaks)] <- "TF"
  sMAF_peaks <- rbind(sMAF_peaks, thesePeaks)
}


cofactor_peaks <- c()
cofactor_peaks <- c()
for (i in 1:nrow(Cofactor_table)) {
  thesePeaks <- read.table(Cofactor_table[i,2], header=F, sep="\t", stringsAsFactors=F)
  thesePeaks <- cbind(thesePeaks, rep(Cofactor_table[i,1], nrow(thesePeaks)))
  colnames(thesePeaks)[ncol(thesePeaks)] <- "TF"
  cofactor_peaks <- rbind(cofactor_peaks, thesePeaks)
}




all_other_TF_peaks <- c()
for (i in 1:nrow(TFs_table)) {
  print(paste(i, nrow(TFs_table)))
  thesePeaks <- read.table(TFs_table[i,2], header=F, sep="\t", stringsAsFactors=F)
  thesePeaks <- cbind(thesePeaks, rep(TFs_table[i,1], nrow(thesePeaks)))
  colnames(thesePeaks)[ncol(thesePeaks)] <- "TF"
  all_other_TF_peaks <- rbind(all_other_TF_peaks, thesePeaks)
}


################################################################################
#For each set, create a dataframe, center the peaks around the summits,
#and extend by 50bp on each side.
################################################################################

sMAF_peaks <- as.data.frame(sMAF_peaks)
colnames(sMAF_peaks)[1:3] <- c("chr", "start", "end")
newStarts <- sMAF_peaks$start + sMAF_peaks$V10 - 50
newEnds <- sMAF_peaks$start + sMAF_peaks$V10 + 50
sMAF_peaks$start <- newStarts
sMAF_peaks$end <- newEnds
sMAF_peaks_gr <- makeGRangesFromDataFrame(sMAF_peaks)


cofactor_peaks <- as.data.frame(cofactor_peaks)
colnames(cofactor_peaks)[1:3] <- c("chr", "start", "end")
newStarts <- cofactor_peaks$start + cofactor_peaks$V10 - 50
newEnds <- cofactor_peaks$start + cofactor_peaks$V10 + 50
cofactor_peaks$start <- newStarts
cofactor_peaks$end <- newEnds
cofactor_peaks_gr <- makeGRangesFromDataFrame(cofactor_peaks)


all_other_TF_peaks <- as.data.frame(all_other_TF_peaks)
colnames(all_other_TF_peaks)[1:3] <- c("chr", "start", "end")
newStarts <- all_other_TF_peaks$start + all_other_TF_peaks$V10 - 50
newEnds <- all_other_TF_peaks$start + all_other_TF_peaks$V10 + 50
all_other_TF_peaks$start <- newStarts
all_other_TF_peaks$end <- newEnds
all_other_TF_peaks_gr <- makeGRangesFromDataFrame(all_other_TF_peaks)


################################################################################
#First, get sMAF only peaks, by overlapping that set with other sets, and removing
#any peak with sufficient overlap (defined as?... 50bp overlap)
################################################################################

smaf_cofactor_overlap <- as.data.frame(findOverlaps(sMAF_peaks_gr, cofactor_peaks_gr, minoverlap=50))
smaf_other_overlap <- as.data.frame(findOverlaps(sMAF_peaks_gr, all_other_TF_peaks_gr, minoverlap=50))

sMAF_only_peaks <- sMAF_peaks[-unique(c(smaf_cofactor_overlap[,1], smaf_other_overlap[,1])),]
sMAF_only_peaks_gr <- makeGRangesFromDataFrame(sMAF_only_peaks)



################################################################################
#Second, for each type of sMAF, determine whether or not it overlaps
#with another sMAF type.  This will define your sMAF only overlapping case.
################################################################################


sMAF_only_peaks_selfOverlap <- as.data.frame(findOverlaps(sMAF_only_peaks_gr, sMAF_only_peaks_gr, minoverlap=50))
sMAF_only_peaks_selfOverlap <- sMAF_only_peaks_selfOverlap[sMAF_only_peaks_selfOverlap[,1]!=sMAF_only_peaks_selfOverlap[,2],]

for (i in 1:nrow(sMAF_only_peaks_selfOverlap)) {
  thisSet <- as.numeric(sMAF_only_peaks_selfOverlap[i,])
  thisSet <- thisSet[order(thisSet)]
  sMAF_only_peaks_selfOverlap[i,] <- thisSet
}
sMAF_only_peaks_selfOverlap <- unique(sMAF_only_peaks_selfOverlap)


sMAF_peaks_overlapping <- data.frame()
for (i in 1:nrow(sMAF_only_peaks_selfOverlap)) {
  firstSet <- sMAF_only_peaks[sMAF_only_peaks_selfOverlap[i,1],]
  secondSet <- sMAF_only_peaks[sMAF_only_peaks_selfOverlap[i,2],]

  theSummit_1 <- mean(as.numeric(firstSet[1,2:3]))
  theSummit_2 <- mean(as.numeric(secondSet[1,2:3]))
  theSummit <- floor(mean(c(theSummit_1, theSummit_2)))

  newStart <- theSummit-50
  newEnd <- theSummit+50

  newSet <- firstSet
  newSet[1,2] <- newStart
  newSet[1,3] <- newEnd
  newSet[1,10] <- 51
  newSet[1,11] <- paste(firstSet[1,11], secondSet[1,11], sep=";")

  sMAF_peaks_overlapping <- rbind(sMAF_peaks_overlapping, newSet)

}

sMAF_peaks_overlapping_gr <- makeGRangesFromDataFrame(sMAF_peaks_overlapping)



overlaps_with_smafsmaf_cCREs <- as.data.frame(findOverlaps(sMAF_peaks_overlapping_gr, cCREs_gr))
the_smaf_overlapping_cCREs <- cCREs[overlaps_with_smafsmaf_cCREs[,2],4]


################################################################################
#Third, grab all of the peaks which are NOT sMAF-sMAF overlapping.
################################################################################

any_sMAF_peak_overlap <- as.data.frame(findOverlaps(sMAF_peaks_gr, sMAF_peaks_gr, minoverlap=50))
any_sMAF_peak_overlap <- any_sMAF_peak_overlap[any_sMAF_peak_overlap[,1]!=any_sMAF_peak_overlap[,2],]


sMAF_peaks_no_selfOverlap <- sMAF_peaks[-unique(any_sMAF_peak_overlap[,1]),]
sMAF_peaks_no_selfOverlap_gr <- makeGRangesFromDataFrame(sMAF_peaks_no_selfOverlap)


################################################################################
#Fourth, get Cofactor overlaps.
################################################################################

sMAF_Cofactor_overlaps <- as.data.frame(findOverlaps(sMAF_peaks_no_selfOverlap_gr, cofactor_peaks_gr, minoverlap=50))

#Make the "peak" list as we did with the smaf self overlaps.


sMAF_Cofactor_peaks <- data.frame()
for (i in 1:nrow(sMAF_Cofactor_overlaps)) {
  firstSet <- sMAF_peaks_no_selfOverlap[sMAF_Cofactor_overlaps[i,1],]
  secondSet <- cofactor_peaks[sMAF_Cofactor_overlaps[i,2],]

  theSummit_1 <- mean(as.numeric(firstSet[1,2:3]))
  theSummit_2 <- mean(as.numeric(secondSet[1,2:3]))
  theSummit <- floor(mean(c(theSummit_1, theSummit_2)))

  newStart <- theSummit-50
  newEnd <- theSummit+50

  newSet <- firstSet
  newSet[1,2] <- newStart
  newSet[1,3] <- newEnd
  newSet[1,10] <- 51
  newSet[1,11] <- paste(firstSet[1,11], secondSet[1,11], sep=";")

  sMAF_Cofactor_peaks <- rbind(sMAF_Cofactor_peaks, newSet)

}

sMAF_Cofactor_peaks_gr <- makeGRangesFromDataFrame(sMAF_Cofactor_peaks)

################################################################################
#Fifth, remove all cases of sMAF overlap, however small, with itself
#or with NEF2.  Then perform overlaps with other TFs, requiring 50bp overlap.
################################################################################

any_sMAF_peak_overlap_strict <- as.data.frame(findOverlaps(sMAF_peaks_gr, sMAF_peaks_gr, minoverlap=1))
any_sMAF_peak_overlap_strict <- any_sMAF_peak_overlap_strict[any_sMAF_peak_overlap_strict[,1]!=any_sMAF_peak_overlap_strict[,2],]


sMAF_peaks_firstFilter <- sMAF_peaks[-unique(any_sMAF_peak_overlap_strict[,1]),]
sMAF_peaks_firstFilter_gr <- makeGRangesFromDataFrame(sMAF_peaks_firstFilter)

any_sMAF_peak_Cofactor_overlap_strict <- as.data.frame(findOverlaps(sMAF_peaks_firstFilter_gr, cofactor_peaks_gr, minoverlap=1))
sMAF_peaks_secondFilter <- sMAF_peaks_firstFilter[-unique(any_sMAF_peak_Cofactor_overlap_strict[,1]),]
sMAF_peaks_secondFilter_gr <- makeGRangesFromDataFrame(sMAF_peaks_secondFilter)

any_sMAF_with_otherPeak_overlap <- as.data.frame(findOverlaps(sMAF_peaks_secondFilter_gr, all_other_TF_peaks_gr, minoverlap=50))

################################################################################
#Make the "peak" list as we did with the smaf self overlaps.
################################################################################
sMAF_Other_peaks <- sMAF_peaks_secondFilter[unique(any_sMAF_with_otherPeak_overlap[,1]),]

sMAF_Other_peaks_gr <- makeGRangesFromDataFrame(sMAF_Other_peaks)


################################################################################
#Sixth, let's get a control set of peaks.
#for the "other" set, remove cases overlapping with sMAFs or cofactors.
#then do a self-overlap.
#make it unique.
#remove cases of the same TF overlapping.
#Sample a number equal to the number of sMAF_Other peaks.
################################################################################

any_smAF_other_overlap <- as.data.frame(findOverlaps(all_other_TF_peaks_gr, sMAF_peaks_gr, minoverlap=1))
all_other_TF_peaks_firstFilter <- all_other_TF_peaks[-unique(any_smAF_other_overlap[,1]),]
all_other_TF_peaks_firstFilter_gr <- makeGRangesFromDataFrame(all_other_TF_peaks_firstFilter)

any_cofactor_overlap <- as.data.frame(findOverlaps(all_other_TF_peaks_firstFilter_gr, cofactor_peaks_gr, minoverlap=1))
all_other_TF_peaks_secondFilter <- all_other_TF_peaks_firstFilter[-unique(any_cofactor_overlap[,1]),]
all_other_TF_peaks_secondFilter_gr <- makeGRangesFromDataFrame(all_other_TF_peaks_secondFilter)

all_other_selfOverlap <- as.data.frame(findOverlaps(all_other_TF_peaks_secondFilter_gr, all_other_TF_peaks_secondFilter_gr, minoverlap=50))
all_other_selfOverlap <- all_other_selfOverlap[all_other_selfOverlap[,1]!=all_other_selfOverlap[,2],]
dim(all_other_selfOverlap[all_other_selfOverlap[,2]>all_other_selfOverlap[,1],])

all_other_selfOverlap <- all_other_selfOverlap[all_other_selfOverlap[,2]>all_other_selfOverlap[,1],]


theSample <- sample(1:nrow(all_other_selfOverlap), nrow(sMAF_Other_peaks), replace=F)

all_other_selfOverlap_sample <- all_other_selfOverlap[theSample,]


all_other_selfOverlap_peaks <- data.frame()
for (i in 1:nrow(all_other_selfOverlap_sample)) {
  firstSet <- all_other_TF_peaks_secondFilter[all_other_selfOverlap_sample[i,1],]
  secondSet <- all_other_TF_peaks_secondFilter[all_other_selfOverlap_sample[i,2],]

  theSummit_1 <- mean(as.numeric(firstSet[1,2:3]))
  theSummit_2 <- mean(as.numeric(secondSet[1,2:3]))
  theSummit <- floor(mean(c(theSummit_1, theSummit_2)))

  newStart <- theSummit-50
  newEnd <- theSummit+50

  newSet <- firstSet
  newSet[1,2] <- newStart
  newSet[1,3] <- newEnd
  newSet[1,10] <- 51
  newSet[1,11] <- paste(firstSet[1,11], secondSet[1,11], sep=";")

  all_other_selfOverlap_peaks <- rbind(all_other_selfOverlap_peaks, newSet)

}

all_other_selfOverlap_peaks <- all_other_selfOverlap_peaks[order(all_other_selfOverlap_peaks[,7], decreasing=T),]



################################################################################
#Save each of these 4 peak sets.
################################################################################

saveFile <- paste(outDir, "sMAF_sMAF_cobound_peaks.bed", sep="")
write.table(sMAF_peaks_overlapping, saveFile, row.names=F, col.names=F, sep="\t", quote=F)


saveFile <- paste(outDir, "sMAF_Cofactor_cobound_peaks.bed", sep="")
write.table(sMAF_Cofactor_peaks, saveFile, row.names=F, col.names=F, sep="\t", quote=F)


saveFile <- paste(outDir, "sMAF_Other_cobound_peaks.bed", sep="")
write.table(sMAF_Other_peaks, saveFile, row.names=F, col.names=F, sep="\t", quote=F)


saveFile <- paste(outDir, "other_other_cobound_peaks.bed", sep="")
write.table(all_other_selfOverlap_peaks, saveFile, row.names=F, col.names=F, sep="\t", quote=F)




################################################################################
#For each of these peak sets, create a GR and then get the expression levels of
#the nearest genes.
################################################################################

sMAF_peaks_overlapping_gr <- makeGRangesFromDataFrame(sMAF_peaks_overlapping)
sMAF_Cofactor_peaks_gr <- makeGRangesFromDataFrame(sMAF_Cofactor_peaks)
sMAF_Other_peaks_gr <- makeGRangesFromDataFrame(sMAF_Other_peaks)
other_other_peaks_gr <- makeGRangesFromDataFrame(all_other_selfOverlap_peaks)

#Gotta get a gene set genomic ranges object first.

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

gene_locs <- data.frame(chr=fullCounts1$Chromosome, start=fullCounts1$TSS, end=fullCounts1$TSS)
gene_locs_gr <- makeGRangesFromDataFrame(gene_locs)



#Get some distance metrics first.

a1 <- as.data.frame(distanceToNearest(sMAF_peaks_overlapping_gr, gene_locs_gr))
a2 <- as.data.frame(distanceToNearest(sMAF_Cofactor_peaks_gr, gene_locs_gr))
a3 <- as.data.frame(distanceToNearest(sMAF_Other_peaks_gr, gene_locs_gr))
a4 <- as.data.frame(distanceToNearest(other_other_peaks_gr, gene_locs_gr))


#Get the expression of nearest, restricted to 20kb away.

nearest_sMAF_sMAF <- as.data.frame(nearest(sMAF_peaks_overlapping_gr, gene_locs_gr))
nearest_sMAF_Cofactor <- as.data.frame(nearest(sMAF_Cofactor_peaks_gr, gene_locs_gr))
nearest_sMAF_Other <- as.data.frame(nearest(sMAF_Other_peaks_gr, gene_locs_gr))
nearest_other_other <- as.data.frame(nearest(other_other_peaks_gr, gene_locs_gr))


nearest_sMAF_sMAF_expr <- unlist(fullCounts1[nearest_sMAF_sMAF[,1],"TPM"])

nearest_sMAF_Cofactor_expr <- unlist(fullCounts1[nearest_sMAF_Cofactor[,1],"TPM"])

nearest_sMAF_Other_expr <- unlist(fullCounts1[nearest_sMAF_Other[,1],"TPM"])

nearest_other_other_expr <- unlist(fullCounts1[nearest_other_other[,1],"TPM"])



#Make some data frames of these.

nearest_sMAF_sMAF_expr <- cbind(a1, nearest_sMAF_sMAF_expr, rep("sMAF_sMAF", length(nearest_sMAF_sMAF_expr)))
colnames(nearest_sMAF_sMAF_expr) <- c("PeakIndex", "GeneIndex", "distance", "TPM", "Type")

nearest_sMAF_Cofactor_expr <- cbind(a2, nearest_sMAF_Cofactor_expr, rep("sMAF_Cofactor", length(nearest_sMAF_Cofactor_expr)))
colnames(nearest_sMAF_Cofactor_expr) <- c("PeakIndex", "GeneIndex", "distance", "TPM", "Type")

nearest_sMAF_Other_expr <- cbind(a3, nearest_sMAF_Other_expr, rep("sMAF_Other", length(nearest_sMAF_Other_expr)))
colnames(nearest_sMAF_Other_expr) <- c("PeakIndex", "GeneIndex", "distance", "TPM", "Type")

nearest_other_other_expr <- cbind(a3, nearest_other_other_expr, rep("Other_Other", length(nearest_other_other_expr)))
colnames(nearest_other_other_expr) <- c("PeakIndex", "GeneIndex", "distance", "TPM", "Type")






#

graphDF <- rbind(nearest_sMAF_sMAF_expr[,c("TPM", "Type")],
  nearest_sMAF_Cofactor_expr[,c("TPM", "Type")],
  nearest_sMAF_Other_expr[,c("TPM", "Type")],
  nearest_other_other_expr[,c("TPM", "Type")])

graphDF$logTPM <- log(graphDF$TPM+0.01)
graphDF$Type <- factor(graphDF$Type, levels=c("sMAF_sMAF", "sMAF_Cofactor", "sMAF_Other", "Other_Other"))



saveFile <- paste(outDir, "Moyers_Figure4D.pdf", sep="")
p <- ggplot(graphDF, aes(x=Type, y=logTPM)) + theme_bw() + geom_boxplot() + ylab("log(TPM) of Nearest Gene") + xlab("") +
  theme(axis.text= element_text(size=15), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5)) + ylim(-5,10) +
  scale_x_discrete(labels=c('sMAF-sMAF', 'sMAF-Cofactor', 'sMAF-Other', 'non-sMAF'))
ggsave(saveFile)




################################################################################
#Now the challenge of deriving motifs for these.
################################################################################




hotSites <- read.delim(hotSites, header=F, sep="\t", stringsAsFactors=F)
hotSites <- as.data.frame(hotSites)
colnames(hotSites)[1:3] <- c("chr", "start", "end")
hotSites_gr <- makeGRangesFromDataFrame(hotSites)



genome <- FaFile(file=genomeFile)




this_outDir <- paste(outDir, "sMAF_sMAF_meme/", sep="")
running_meme_pipeline(thePeaks=sMAF_peaks_overlapping, hotSites_gr, genome, this_outDir)


this_outDir <- paste(outDir, "sMAF_Cofactor_meme/", sep="")
running_meme_pipeline(thePeaks=sMAF_Cofactor_peaks, hotSites_gr, genome, this_outDir)

this_outDir <- paste(outDir, "sMAF_Other_meme/", sep="")
running_meme_pipeline(thePeaks=sMAF_Other_peaks, hotSites_gr, genome, this_outDir)


this_outDir <- paste(outDir, "Other_Other_meme/", sep="")
running_meme_pipeline(thePeaks=all_other_selfOverlap_peaks, hotSites_gr, genome, this_outDir)





################################################################################
#Intersect with 3d interaction data to identify support of connections for each set.
################################################################################



hiccup <- "/cluster/home/bmoyers/ENCODE_500Plus/Hiccups_3dContact_data_2023April04/ENCFF264RQT.bedpe.gz"
hiccup <- read.table(hiccup, header=F, sep="\t", stringsAsFactors=F)

hiccup_first <- as.data.frame(hiccup[,1:3])
colnames(hiccup_first) <- c("chr", "start", "end")
hiccup_first_gr <- makeGRangesFromDataFrame(hiccup_first)
a <- width(hiccup_first_gr)
hiccup_first_gr <- hiccup_first_gr[which(a==1001)]
hiccup_highRes_1 <- hiccup[which(a==1001),]

hiccup_second <- as.data.frame(hiccup[,4:6])
colnames(hiccup_second) <- c("chr", "start", "end")
hiccup_second_gr <- makeGRangesFromDataFrame(hiccup_second)
a <- width(hiccup_second_gr)
hiccup_second_gr <- hiccup_second_gr[which(a==1001)]
hiccup_highRes_2 <- hiccup[which(a==1001),]


hiccup_first <- hiccup_highRes_1[,c(1:3)]
colnames(hiccup_first) <- c("chr", "start", "end")

hiccup_second <- as.data.frame(hiccup_highRes_2[,4:6])
colnames(hiccup_second) <- c("chr", "start", "end")

hiccup_highRes_both <- rbind(hiccup_first, hiccup_second)
hiccup_highRes_both <- unique(hiccup_highRes_both)
hiccup_highRes_both_gr <- makeGRangesFromDataFrame(hiccup_highRes_both)


################################################################################
#Grab overlap with each of the 4 sets.
################################################################################

intersect_hiccup_sMAFsMAF <- as.data.frame(findOverlaps(sMAF_peaks_overlapping_gr, hiccup_highRes_both_gr))

intersect_hiccup_sMAFCofactor <- as.data.frame(findOverlaps(sMAF_Cofactor_peaks_gr, hiccup_highRes_both_gr))

intersect_hiccup_sMAFOther <- as.data.frame(findOverlaps(sMAF_Other_peaks_gr, hiccup_highRes_both_gr))

intersect_hiccup_OtherOther <- as.data.frame(findOverlaps(all_other_selfOverlap_peaks_gr, hiccup_highRes_both_gr))


hiccup_overlap_data <- rbind(c(length(unique(intersect_hiccup_sMAFsMAF[,1])), length(sMAF_peaks_overlapping_gr)),
      c(length(unique(intersect_hiccup_sMAFCofactor[,1])), length(sMAF_Cofactor_peaks_gr)),
      c(length(unique(intersect_hiccup_sMAFOther[,1])), length(sMAF_Other_peaks_gr)),
      c(length(unique(intersect_hiccup_OtherOther[,1])), length(all_other_selfOverlap_peaks_gr)))

hiccup_overlap_data <- as.data.frame(hiccup_overlap_data)
colnames(hiccup_overlap_data) <- c("regions_in_hiccups", "total_regions")
rownames(hiccup_overlap_data) <- c("sMAFsMAF", "sMAFCofactor", "sMAFOther", "OtherOther")
hiccup_overlap_data$Fraction <- hiccup_overlap_data$regions_in_hiccups / hiccup_overlap_data$total_regions

hiccup_overlap_data$set <- rownames(hiccup_overlap_data)
hiccup_overlap_data$set <- factor(hiccup_overlap_data$set, levels=hiccup_overlap_data$set)


################################################################################
#Make a barpot of these sets.
################################################################################

saveFile <- paste(outDir, "Moyers_Figure4E.pdf", sep="")
p <- ggplot(data=hiccup_overlap_data, aes(x=set, y=Fraction)) + theme_bw() +
  geom_bar(stat="identity") +  theme(axis.text= element_text(size=15), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5)) +
  ylab("Fraction with Connection") + xlab("") +
  scale_x_discrete(labels=c('sMAF-sMAF', 'sMAF-Cofactor', 'sMAF-Other', 'non-sMAF'))
ggsave(saveFile)




#
