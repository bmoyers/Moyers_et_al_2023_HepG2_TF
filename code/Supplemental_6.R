#!/usr/bin/env R
#Supplemental_6.R

################################################################################
#This script takes a nested directory of gkmpredict scores for various regions
#using gkmsvm models for all ChIP-seq datasets to produce supplemental figure 6.
#These pre-computed scores can be generated using the script Running_gkmpredict.R
#
#Run under R 4.1.0
#This script takes as arguments:
# outDir: Path where files are to be written.  Note that it is expected that this
#     outDir argument is the same as was used in Running_gkmpredict.R
# unbound_Bed: cCRE v4 regions which are in open ATAC-seq peaks in HepG2
#     and which did not contain a ChIP-seq peak in any of our datasets.  Provided
#     as unbound_regions_cCREs_in_ATAC.bed
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
#
################################################################################

################################################################################
################################################################################
#Load Libraries.
################################################################################
################################################################################

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
exprDir <- args[3]


theExprs <- list.files(exprDir, full.names=T, pattern="Preferred")
theExprs <- theExprs[-grep("H3K[3|4]", theExprs)]

gkm_outDir <- paste(outDir, "gkmsvm_outputs/", sep="")

outDir_Bound <- paste(gkm_outDir, "gkm_scores_Bound_v4cCREs/", sep="")
outDir_Unbound_cCREs <- paste(gkm_outDir, "gkm_scores_Unbound_cCREs_v4cCREs/", sep="")
outDir_Matched <- paste(gkm_outDir, "gkm_scores_Matched_v4cCREs/", sep="")
outDir_cCREs_Bound <- paste(gkm_outDir, "gkm_scores_cCREs_Bound_v4cCREs/", sep="")
outDir_cCREs_Unbound <- paste(gkm_outDir, "gkm_scores_cCREs_Unbound_v4cCREs/", sep="")



cCREs_unbound_theMax <- c()
matched_theMax <- c()
all_cCRE_bound <- c()

for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisName <- strsplit(theExprs[i], split="/")[[1]][length(strsplit(theExprs[i], split="/")[[1]])]
  thisName <- gsub(".bed.gz", "", thisName)
  thisBed <- theExprs[i]

  this_cCRE_bound <- paste(outDir_cCREs_Bound, thisName, ".txt", sep="")
  this_Matched <- paste(outDir_Matched, thisName, ".txt", sep="")
  this_Unbound <- paste(outDir_Unbound_cCREs, thisName, ".txt", sep="")
  this_cCRE_unbound <- paste(outDir_cCREs_Unbound, thisName, ".txt", sep="")

  if(file.exists(this_cCRE_bound) && file.exists(this_Matched) && file.exists(this_Unbound) && file.exists(this_cCRE_unbound)) {

    this_cCRE_bound <- read.delim(this_cCRE_bound, header=F, sep="\t", stringsAsFactors=F)
    this_cCRE_bound_sd <- sd(this_cCRE_bound[,2])
    this_cCRE_bound_mean <- mean(this_cCRE_bound[,2])

    this_cCRE_bound <- as.data.frame(this_cCRE_bound)
    this_cCRE_bound$Score_BoundcCREsNormalized <- (this_cCRE_bound[,2]-this_cCRE_bound_mean)/this_cCRE_bound_sd


    this_Matched <- read.delim(this_Matched, header=F, sep="\t", stringsAsFactors=F)
    this_Matched <- as.data.frame(this_Matched)
    this_Matched$Score_BoundcCREsNormalized <- (this_Matched[,2]-this_cCRE_bound_mean)/this_cCRE_bound_sd

    this_Unbound <- read.delim(this_Unbound, header=F, sep="\t", stringsAsFactors=F)
    this_Unbound <- as.data.frame(this_Unbound)
    this_Unbound$Category <- rep("Unbound_cCREs", nrow(this_Unbound))
    this_Unbound$Score_BoundcCREsNormalized <- (this_Unbound[,2]-this_cCRE_bound_mean)/this_cCRE_bound_sd

    if(i==1) {
      cCREs_unbound_theMax <- c(this_Unbound$Score_BoundcCREsNormalized)
      matched_theMax <- c(this_Matched$Score_BoundcCREsNormalized)
      all_cCRE_bound <- c(all_cCRE_bound, this_cCRE_bound$Score_BoundcCREsNormalized)
    }


    if(i>1) {
      cCREs_unbound_theMax <- rowMaxs(cbind(cCREs_unbound_theMax, this_Unbound$Score_BoundcCREsNormalized))
      matched_theMax <- rowMaxs(cbind(matched_theMax, this_Matched$Score_BoundcCREsNormalized))
      all_cCRE_bound <- c(all_cCRE_bound, this_cCRE_bound$Score_BoundcCREsNormalized)
    }

  }

}



all_cCRE_bound <- cbind(all_cCRE_bound, rep("Bound_cCREs_w_TF", length(all_cCRE_bound)))
matched_theMax <- cbind(matched_theMax, rep("Matched_Null_MaximumScore", length(matched_theMax)))
cCREs_unbound_theMax <- cbind(cCREs_unbound_theMax, rep("Unbound_cCREs_MaximumScore", length(cCREs_unbound_theMax)))

colnames(all_cCRE_bound) <- c("Score", "Category")
colnames(matched_theMax) <- c("Score", "Category")
colnames(cCREs_unbound_theMax) <- c("Score", "Category")

graphDF <- rbind(all_cCRE_bound, matched_theMax, cCREs_unbound_theMax)
graphDF <- as.data.frame(graphDF)
graphDF$Score <- as.numeric(as.character(graphDF$Score))
graphDF$Category <- factor(graphDF$Category, levels=unique(graphDF$Category))


saveFile <- paste(outDir, "Max_Scores_matrix.txt", sep="")
write.table(graphDF, saveFile, row.names=F, col.names=T, sep="\t", quote=F)


saveFile <- paste(outDir, "Supplemental_6.pdf", sep="")
p <-  ggplot(graphDF, aes(x=Category, y=Score)) + geom_violin() +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("") + ylab("Normalized gkmsvm score") + xlab("") +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20))
ggsave(saveFile)
