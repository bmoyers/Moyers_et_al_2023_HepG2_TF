#!/usr/bin/env R
#Running_gkmpredict.R

################################################################################
#This script takes ChIP-seq datasets and runs gkmpredict for a gkmsvm model of
#each experiment to calculate scores for various genomic regions sets.
#This is required for the construction of supplemental figure 6.
#
#Note that this script expects the following to be available in the User's PATH:
# python version 2.7.15
# bedtools
# lsgkm (see: https://github.com/Dongwon-Lee/lsgkm)
# nullseq_generate.py script from kmer-SVM (see: https://beerlab.org/kmersvm/)
#This script also constructs several commands which expect an slurm job scheduler
#to be available.  Note that if your local compute environment does not have
#such a scheduler, you will need to rewrite those commands.  Note that
#running these jobs in sequence for each experiment without the use of a job
#scheduler will take an extremely long time.
#
#Run under R 4.1.0
#This script takes as arguments:
# outDir: Path where files are to be written
# unbound_Bed: cCRE v4 regions which are in open ATAC-seq peaks in HepG2
#     and which did not contain a ChIP-seq peak in any of our datasets.  Provided
#     as unbound_regions_cCREs_in_ATAC.bed
# exclusionList: Path to the ENCODE exclusion list, under accession ENCFF356LFX
# cCREs: Path to the Version 4 cCREs provided by Jill Moore in June 2022.
# ATAC_seq: Path to ATAC-seq peaks, under accession ENCFF439EIO
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# genomeFile: hg38 genome fasta file.
# nullseq_indices: nullseq indices for hg38
# gkmsvm_models: a folder containing each of the pre-generated gkmsvm models,
#     provided in compressed format as gkmsvm_models.tar.gz Note that this must be
#     uncompressed for use before this script.  Provided models were constructed
#     using the lsgkm_execution.sh script.
#
################################################################################

################################################################################
################################################################################
#Load Libraries.
################################################################################
################################################################################

library(GenomicRanges)
library(Rsamtools)
library(Biostrings)
library(matrixStats)
library(ggplot2)

################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
unbound_Bed <- args[2]
exclusionList <- args[3]
cCREs <- args[4]
ATAC_seq <- args[5]
exprDir <- args[6]
genomeFile <- args[7]
nullseq_indices <- args[8]
gkmsvm_models <- args[9]


################################################################################
#Go ahead and determine what % of the regions are bound by a peak in non-preferred.
################################################################################


unbound_Bed <- read.delim(unbound_Bed, header=T, sep="\t", stringsAsFactors=F)
unbound_Bed <- unbound_Bed[order(as.character(unbound_Bed[,1]), as.numeric(as.character(unbound_Bed[,2]))),]
unbound_Bed <- as.data.frame(unbound_Bed)
colnames(unbound_Bed)[1:3] <- c("chr", "start", "end")
unbound_Bed_gr <- makeGRangesFromDataFrame(unbound_Bed, ignore.strand=T, keep.extra.columns=F)
unbound_Bed_gr <- keepStandardChromosomes(unbound_Bed_gr, pruning.mode="coarse")



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


cCREs_Intersect <- as.data.frame(findOverlaps(cCREs_gr, exclusionList_gr))
if(nrow(cCREs_Intersect)>0) {
  cCREs <- cCREs[-unique(cCREs_Intersect[,1]),]
  cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=T, keep.extra.columns=F)
}


ATAC_Intersect <- as.data.frame(findOverlaps(cCREs_gr, ATAC_seq_gr))
if(nrow(ATAC_Intersect)>0) {
  cCREs <- cCREs[unique(ATAC_Intersect[,1]),]
  cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=T, keep.extra.columns=F)
}

bound_intersect <- as.data.frame(findOverlaps(cCREs_gr, unbound_Bed_gr))
bound_cCREs <- cCREs
bound_cCREs_gr <- cCREs_gr
if(nrow(cCREs_Intersect)>0) {
  bound_cCREs <- bound_cCREs[-unique(cCREs_Intersect[,1]),]
  bound_cCREs_gr <- makeGRangesFromDataFrame(bound_cCREs, ignore.strand=T, keep.extra.columns=F)
}



theExprs <- list.files(exprDir, full.names=T, pattern="Non-preferred")
theExprs <- theExprs[-grep("H3K[3|4]", theExprs)]

numBound <- rep(0, nrow(unbound_Bed))

for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisFile <- read.delim(theExprs[i], header=F, sep="\t", stringsAsFactors=F)
  colnames(thisFile)[1:3] <- c("chr", "start", "end")
  thisFile_gr <- makeGRangesFromDataFrame(thisFile, ignore.strand=T)

  theIntersect <- as.data.frame(findOverlaps(unbound_Bed_gr, thisFile_gr))

  numBound[unique(theIntersect[,1])] <- numBound[unique(theIntersect[,1])]+1
}

length(numBound)
length(numBound[numBound>0])
length(numBound[numBound>0])/length(numBound)





################################################################################
#Define Genome File.
################################################################################

genome <- FaFile(file=genomeFile)



################################################################################
#Save a FASTA sequence of the unbound cCREs in ATAC peaks
################################################################################
unbound_Bed_names <- paste(unbound_Bed[,1], ":", unbound_Bed[,2], "-", unbound_Bed[,3], sep="")
theSequences <- getSeq(genome, unbound_Bed_gr, as.character=T)
theSequences <- toupper(as.character(theSequences))
dna_theSequences <- DNAStringSet(theSequences)
names(dna_theSequences) <- unbound_Bed_names
saveFile_1 <- paste(outDir, "Unbound_cCRE_sequences.fasta", sep="")
writeXStringSet(dna_theSequences, saveFile_1, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")



################################################################################
#Identify background peak set, remove any overlapping with cCREs.
################################################################################


positive_set <- paste(outDir, "Unbound_cCREs.bed", sep="")
write.table(unbound_Bed, positive_set, row.names=F, col.names=F, sep="\t", quote=F)
negative_set <- paste(outDir, "Matched_regions_for_unbound_cCREs.bed", sep="")

theCommand <- paste("python nullseq_generate.py -G -x 5 -m 1000 -r 1 -o", negative_set, " ", positive_set,  " hg38 ", nullseq_indices, sep="")
system(theCommand)

matched_Bed <- read.delim(negative_set, header=F, sep="\t", stringsAsFactors=F)
matched_Bed <- as.data.frame(matched_Bed)
colnames(matched_Bed) <- c("chr", "start", "end")
matched_Bed_gr <- makeGRangesFromDataFrame(matched_Bed, ignore.strand=T, keep.extra.columns=F)


theIntersect <- as.data.frame(findOverlaps(matched_Bed_gr, cCREs_gr))
if(nrow(theIntersect)>0) {
  matched_Bed <- matched_Bed[-unique(theIntersect[,1]),]
  matched_Bed_gr <- makeGRangesFromDataFrame(matched_Bed, ignore.strand=T, keep.extra.columns=F)
}

summary(unbound_Bed[,3]-unbound_Bed[,2]+1)
summary(matched_Bed[,3]-matched_Bed[,2]+1)

matched_Bed_gr <- keepStandardChromosomes(matched_Bed_gr, pruning.mode="coarse")
matched_Bed <- as.data.frame(matched_Bed_gr)

matched_Bed_names <- paste(matched_Bed[,1], ":", matched_Bed[,2], "-", matched_Bed[,3], sep="")
theSequences <- getSeq(genome, matched_Bed_gr, as.character=T)
theSequences <- toupper(as.character(theSequences))
dna_theSequences <- DNAStringSet(theSequences)
names(dna_theSequences) <- matched_Bed_names
saveFile_1 <- paste(outDir, "Matched_sequences.fasta", sep="")
writeXStringSet(dna_theSequences, saveFile_1, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")





################################################################################
#Now, run gkmpredict on these regions with every single experiment we've got a
#model for.
################################################################################

theExprs <- list.files(exprDir, full.names=T, pattern="Preferred")
theExprs <- theExprs[-grep("H3K[3|4]", theExprs)]


################################################################################
#Write some code to replace the memeDirs situation. Make a gkm_outDir within the
#outDir.  For each experiment, make a directory with the experiment's name.
#Use that location to write out each of the models.
################################################################################

gkm_outDir <- paste(outDir, "gkmsvm_outputs/", sep="")
system(paste("mkdir ", gkm_outDir, sep=""))


fastas_Bound <- paste(gkm_outDir, "gkm_fastas_Bound_v4cCREs/", sep="")
system(paste("mkdir ", fastas_Bound, sep=""))
outDir_Bound <- paste(gkm_outDir, "gkm_scores_Bound_v4cCREs/", sep="")
system(paste("mkdir ", outDir_Bound, sep=""))
outDir_Unbound_cCREs <- paste(gkm_outDir, "gkm_scores_Unbound_cCREs_v4cCREs/", sep="")
system(paste("mkdir ", outDir_Unbound_cCREs, sep=""))
outDir_Matched <- paste(gkm_outDir, "gkm_scores_Matched_v4cCREs/", sep="")
system(paste("mkdir ", outDir_Matched, sep=""))
outDir_cCREs_Bound <- paste(gkm_outDir, "gkm_scores_cCREs_Bound_v4cCREs/", sep="")
system(paste("mkdir ", outDir_cCREs_Bound, sep=""))
outDir_cCREs_Unbound <- paste(gkm_outDir, "gkm_scores_cCREs_Unbound_v4cCREs/", sep="")
system(paste("mkdir ", outDir_cCREs_Unbound, sep=""))


testSeq_Unbound_cCREs <- paste(gkm_outDir, "Unbound_cCRE_sequences.fasta", sep="")
system(paste("mkdir ", testSeq_Unbound_cCREs, sep=""))
testSeq_Matched <- paste(gkm_outDir, "Matched_sequences.fasta", sep="")
system(paste("mkdir ", testSeq_Matched, sep=""))

logFile_bound <- paste(outDir, "gkm_logFile_bound.txt", sep="")
logFile_unbound <- paste(outDir, "gkm_logFile_unbound.txt", sep="")
logFile_matched <- paste(outDir, "gkm_logFile_matched.txt", sep="")
logFile_cCRE_bound <- paste(outDir, "gkm_logFile_cCRE_bound.txt", sep="")
logFile_cCRE_unbound <- paste(outDir, "gkm_logFile_cCRE_unbound.txt", sep="")

for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisName <- strsplit(theExprs[i], split="/")[[1]][length(strsplit(theExprs[i], split="/")[[1]])]
  thisName <- gsub(".bed.gz", "", thisName)
  thisBed <- theExprs[i]
  thisModel <- paste()
  thisDir <- memeDirs[grep(thisName, memeDirs)]
  thisModel <- paste(gkmsvm_models, "/", thisName, ".txt", sep="")

  logFile_bound <- paste(outDir, "gkm_logfiles/", thisName, "_bound.txt", sep="")
  logFile_unbound <- paste(outDir, "gkm_logfiles/", thisName, "_unbound.txt", sep="")
  logFile_matched <- paste(outDir, "gkm_logfiles/", thisName, "_matched.txt", sep="")
  logFile_cCRE_bound <- paste(outDir, "gkm_logfiles/", thisName, "_cCRE_bound.txt", sep="")
  logFile_cCRE_unbound <- paste(outDir, "gkm_logfiles/", thisName, "_cCRE_unbound.txt", sep="")

  if(file.exists(thisModel)) {
    testSeq_Bound <- paste(fastas_Bound, thisName, ".fasta", sep="")
    thisBed <- read.delim(thisBed, header=F, sep="\t", stringsAsFactors=F)
    thisBed <- as.data.frame(thisBed)
    colnames(thisBed) <- c("chr", "start", "end")
    thisBed <- thisBed[order(thisBed[,1], thisBed[,2]),]
    thisBed_gr <- makeGRangesFromDataFrame(thisBed, ignore.strand=T, keep.extra.columns=F)
    thisBed_names <- paste(thisBed[,1], ":", thisBed[,2], "-", thisBed[,3], sep="")
    theSequences <- getSeq(genome, thisBed_gr, as.character=T)
    theSequences <- toupper(as.character(theSequences))
    dna_theSequences <- DNAStringSet(theSequences)
    names(dna_theSequences) <- thisBed_names
    writeXStringSet(dna_theSequences, testSeq_Bound, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")


    thisTF_intersection <- as.data.frame(findOverlaps(bound_cCREs_gr, thisBed_gr))
    if(nrow(thisTF_intersection)>0) {
      this_TF_cCREs_bound <- bound_cCREs[unique(thisTF_intersection[,1]),]
      this_TF_cCREs_bound_gr <- makeGRangesFromDataFrame(this_TF_cCREs_bound, ignore.strand=T, keep.extra.columns=F)
      this_TF_cCREs_unbound <- bound_cCREs[-unique(thisTF_intersection[,1]),]
      this_TF_cCREs_unbound_gr <- makeGRangesFromDataFrame(this_TF_cCREs_unbound, ignore.strand=T, keep.extra.columns=F)

      this_TF_cCREs_bound_fasta <- paste(fastas_Bound, thisName, "_cCREs_bound.fasta", sep="")
      this_TF_cCREs_bound_names <- paste(this_TF_cCREs_bound[,1], ":", this_TF_cCREs_bound[,2], "-", this_TF_cCREs_bound[,3], sep="")
      theSequences <- getSeq(genome, this_TF_cCREs_bound_gr, as.character=T)
      theSequences <- toupper(as.character(theSequences))
      dna_theSequences <- DNAStringSet(theSequences)
      names(dna_theSequences) <- this_TF_cCREs_bound_names
      writeXStringSet(dna_theSequences, this_TF_cCREs_bound_fasta, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")


      this_TF_cCREs_unbound_fasta <- paste(fastas_Bound, thisName, "_cCREs_unbound.fasta", sep="")
      this_TF_cCREs_unbound_names <- paste(this_TF_cCREs_unbound[,1], ":", this_TF_cCREs_unbound[,2], "-", this_TF_cCREs_unbound[,3], sep="")
      theSequences <- getSeq(genome, this_TF_cCREs_unbound_gr, as.character=T)
      theSequences <- toupper(as.character(theSequences))
      dna_theSequences <- DNAStringSet(theSequences)
      names(dna_theSequences) <- this_TF_cCREs_unbound_names
      writeXStringSet(dna_theSequences, this_TF_cCREs_unbound_fasta, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")


      outFile_Bound <- paste(outDir_Bound, thisName, ".txt", sep="")
      firstCommand <- paste("module load g/lsgkm; gkmpredict -T 16", testSeq_Bound, thisModel, outFile_Bound, sep=" ")
      firstCommand <- paste("sbatch --mem=40000 --time=10:00:00 --output=", logFile_bound, " --job-name=gkmsvm --wrap=\"", firstCommand, "\"", sep="")
      system(firstCommand)

      outFile_Unbound_cCREs <- paste(outDir_Unbound_cCREs, thisName, ".txt", sep="")
      secondCommand <- paste("module load g/lsgkm; gkmpredict -T 16", testSeq_Unbound_cCREs, thisModel, outFile_Unbound_cCREs, sep=" ")
      secondCommand <- paste("sbatch --mem=40000 --time=10:00:00 --output=", logFile_unbound, " --job-name=gkmsvm --wrap=\"", secondCommand, "\"", sep="")
      system(secondCommand)

      outFile_Matched <- paste(outDir_Matched, thisName, ".txt", sep="")
      thirdCommand <- paste("module load g/lsgkm; gkmpredict -T 16", testSeq_Matched, thisModel, outFile_Matched, sep=" ")
      thirdCommand <- paste("sbatch --mem=40000 --time=10:00:00 --output=", logFile_matched, " --job-name=gkmsvm --wrap=\"", thirdCommand, "\"", sep="")
      system(thirdCommand)

      outFile_cCREs_Bound <- paste(outDir_cCREs_Bound, thisName, ".txt", sep="")
      fourthCommand <- paste("module load g/lsgkm; gkmpredict -T 16", this_TF_cCREs_bound_fasta, thisModel, outFile_cCREs_Bound, sep=" ")
      fourthCommand <- paste("sbatch --mem=40000 --time=10:00:00 --output=", logFile_cCRE_bound, " --job-name=gkmsvm --wrap=\"", fourthCommand, "\"", sep="")
      system(fourthCommand)

      outFile_cCREs_Unound <- paste(outDir_cCREs_Unbound, thisName, ".txt", sep="")
      fifthCommand <- paste("module load g/lsgkm; gkmpredict -T 16", this_TF_cCREs_unbound_fasta, thisModel, outFile_cCREs_Unound, sep=" ")
      fifthCommand <- paste("sbatch --mem=40000 --time=10:00:00 --output=", logFile_cCRE_unbound, " --job-name=gkmsvm --wrap=\"", fifthCommand, "\"", sep="")
      system(fifthCommand)
    }

  }
}


#
