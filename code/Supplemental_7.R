#!/usr/bin/env R
#Supplemental_7.R

################################################################################
#This script takes a nested directory of gkmpredict scores for various regions
#using gkmsvm models for all ChIP-seq datasets to produce supplemental figure 7.
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
# cCREs: Path to the Version 4 cCREs provided by Jill Moore in June 2022.
# bigWig_path: path to a directory containing bigWigs for select experiments.
#
################################################################################

################################################################################
################################################################################
#Load Libraries.
################################################################################
################################################################################

library(matrixStats)
library(ggplot2)

options("scipen"=100)

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

#gkm_outDir <- paste(outDir, "gkmsvm_outputs/", sep="")
gkm_outDir <- outDir

outDir_Bound <- paste(gkm_outDir, "gkm_scores_Bound_v4cCREs/", sep="")
outDir_Unbound_cCREs <- paste(gkm_outDir, "gkm_scores_Unbound_cCREs_v4cCREs/", sep="")
outDir_Matched <- paste(gkm_outDir, "gkm_scores_Matched_v4cCREs/", sep="")
outDir_cCREs_Bound <- paste(gkm_outDir, "gkm_scores_cCREs_Bound_v4cCREs/", sep="")
outDir_cCREs_Unbound <- paste(gkm_outDir, "gkm_scores_cCREs_Unbound_v4cCREs/", sep="")

cCREs <- read.table(cCREs, header=F, sep="\t", stringsAsFactors=F)
cCREs$name <- paste(cCREs[,1], ":", cCREs[,2], "-", cCREs[,3], sep="")


cCREs_unbound_theMax <- c()
matched_theMax <- c()
all_cCRE_bound <- c()
max_cCRE_bound <- rep(NA, nrow(cCREs))
matched_all <- c()
cCREs_unbound_all <- c()

for (i in 1:length(theExprs)) {
  print(paste(i, length(theExprs)))
  thisName <- strsplit(theExprs[i], split="/")[[1]][length(strsplit(theExprs[i], split="/")[[1]])]
  thisName <- gsub(".bed.gz", "", thisName)
  thisBed <- theExprs[i]

  this_cCRE_bound <- paste(outDir_cCREs_Bound, thisName, ".txt", sep="")
  this_Matched <- paste(outDir_Matched, thisName, ".txt", sep="")
  this_Unbound <- paste(outDir_Unbound_cCREs, thisName, ".txt", sep="")


  if(file.exists(this_cCRE_bound) && file.exists(this_Matched) && file.exists(this_Unbound)) {

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
      matched_all <- c(matched_all, this_Matched$Score_BoundcCREsNormalized)
      cCREs_unbound_all <- c(cCREs_unbound_all, this_Unbound$Score_BoundcCREsNormalized)


      toCheck <- match(this_cCRE_bound[,1], cCREs$name)
      toCheck_nas <- toCheck[is.na(max_cCRE_bound[toCheck])]
      toCheck_num <- toCheck[!is.na(max_cCRE_bound[toCheck])]
      max_cCRE_bound[toCheck_nas] <- this_cCRE_bound[which(toCheck%in%toCheck_nas),3]
      max_cCRE_bound[toCheck_num] <- rowMaxs(cbind(this_cCRE_bound[which(toCheck%in%toCheck_num),3], max_cCRE_bound[toCheck_num]))

    }


    if(i>1) {
      cCREs_unbound_theMax <- rowMaxs(cbind(cCREs_unbound_theMax, this_Unbound$Score_BoundcCREsNormalized))
      matched_theMax <- rowMaxs(cbind(matched_theMax, this_Matched$Score_BoundcCREsNormalized))
      all_cCRE_bound <- c(all_cCRE_bound, this_cCRE_bound$Score_BoundcCREsNormalized)
      matched_all <- c(matched_all, this_Matched$Score_BoundcCREsNormalized)
      cCREs_unbound_all <- c(cCREs_unbound_all, this_Unbound$Score_BoundcCREsNormalized)


      toCheck <- match(this_cCRE_bound[,1], cCREs$name)
      toCheck_nas <- toCheck[is.na(max_cCRE_bound[toCheck])]
      toCheck_num <- toCheck[!is.na(max_cCRE_bound[toCheck])]
      max_cCRE_bound[toCheck_nas] <- this_cCRE_bound[which(toCheck%in%toCheck_nas),3]
      max_cCRE_bound[toCheck_num] <- rowMaxs(cbind(this_cCRE_bound[which(toCheck%in%toCheck_num),3], max_cCRE_bound[toCheck_num]))

    }

  }

}



all_cCRE_bound <- cbind(all_cCRE_bound, rep("Bound cCREs, All Scores for bound DAPs", length(all_cCRE_bound)))
matched_theMax <- cbind(matched_theMax, rep("Nullmatched sequences, Max Scores for all DAPs", length(matched_theMax)))
cCREs_unbound_theMax <- cbind(cCREs_unbound_theMax, rep("Unbound cCREs, Max Scores for all DAPs", length(cCREs_unbound_theMax)))
matched_all <- cbind(matched_all, rep("Nullmatched sequences, All Scores for all DAPs", length(matched_all)))
cCREs_unbound_all <- cbind(cCREs_unbound_all, rep("Unbound cCREs, All Scores for all DAPs", length(cCREs_unbound_all)))
max_cCRE_bound <- cbind(max_cCRE_bound, rep("Bound cCREs, Max Scores for bound DAPs", length(max_cCRE_bound)))
max_cCRE_bound <- max_cCRE_bound[!is.na(max_cCRE_bound[,1]),]

colnames(all_cCRE_bound) <- c("Score", "Category")
colnames(matched_theMax) <- c("Score", "Category")
colnames(cCREs_unbound_theMax) <- c("Score", "Category")
colnames(matched_all) <- c("Score", "Category")
colnames(cCREs_unbound_all) <- c("Score", "Category")
colnames(max_cCRE_bound) <- c("Score", "Category")

graphDF <- rbind(all_cCRE_bound, max_cCRE_bound, matched_all, matched_theMax, cCREs_unbound_all, cCREs_unbound_theMax)
graphDF <- as.data.frame(graphDF)
graphDF$Score <- as.numeric(as.character(graphDF$Score))
graphDF$Category <- factor(graphDF$Category, levels=c("Bound cCREs, All Scores for bound DAPs",
  "Nullmatched sequences, All Scores for all DAPs", "Unbound cCREs, All Scores for all DAPs",
  "Bound cCREs, Max Scores for bound DAPs", "Nullmatched sequences, Max Scores for all DAPs",
  "Unbound cCREs, Max Scores for all DAPs"))


outDir <- "/cluster/home/bmoyers/Figures/ENCODE_500plus/Results_post_July5/PaperFigures_2023August22/"
saveFile <- paste(outDir, "GKMSVM_Scores_matrix.txt", sep="")
write.table(graphDF, saveFile, row.names=F, col.names=T, sep="\t", quote=F)


saveFile <- paste(outDir, "Supplemental_7.pdf", sep="")
p <-  ggplot(graphDF, aes(x=Category, y=Score)) + geom_violin() + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("") + ylab("Normalized gkmsvm score") + xlab("") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=15))
ggsave(saveFile)



thecategories <- c("Bound cCREs, All Scores for bound DAPs",
  "Nullmatched sequences, All Scores for all DAPs", "Unbound cCREs, All Scores for all DAPs",
  "Bound cCREs, Max Scores for bound DAPs", "Nullmatched sequences, Max Scores for all DAPs",
  "Unbound cCREs, Max Scores for all DAPs")

pvalTable <- matrix(nrow=length(thecategories), ncol=length(thecategories))
for (i in 1:(length(thecategories)-1)) {
  firstSet <- graphDF[graphDF$Category==thecategories[i],]
  for (j in (i+1):length(thecategories)) {
    secondSet <- graphDF[graphDF$Category==thecategories[j],]
    thisTest <- t.test(firstSet$Score, secondSet$Score, alternative="two.sided")
    pvalTable[i,j] <- thisTest$p.value
  }
}

colnames(pvalTable) <- thecategories
rownames(pvalTable) <- thecategories

saveFile <- paste(outDir, "Supplemental_7_table.txt", sep="")
write.table(pvalTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)


meanTable <- c()
for (i in 1:length(thecategories)) {
  firstSet <- graphDF[graphDF$Category==thecategories[i],]
  meanTable <- rbind(meanTable, c(thecategories[i], mean(as.numeric(firstSet$Score))))
}


meanTable[,1] <- as.character(thecategories)
colnames(meanTable) <- c("Category", "Mean Normalized Score")

saveFile <- paste(outDir, "Supplemental_7_table_means.txt", sep="")
write.table(meanTable, saveFile, row.names=T, col.names=T, sep="\t", quote=F)


################################################################################
#Testing a subset of bigWigs for curious reviewers.
################################################################################

#cd /cluster/home/bmoyers/ENCODE_500Plus/deepTools

#source ~/miniconda3/bin/activate
#conda activate deepTools

outDir <- "/cluster/home/bmoyers/Figures/ENCODE_500plus/Results_post_July5/PaperFigures_2023August22/"

temp_outDir <- paste(outDir, "temp_deepTools_output", sep="")
if(!file.exists(temp_outDir)) { system(paste("mkdir ", temp_outDir, sep=""))}

figure_outDir <- paste(outDir, "deepTools_falseNegatives_output", sep="")
if(!file.exists(figure_outDir)) { system(paste("mkdir ", figure_outDir, sep=""))}


bigWig_files <- list.files(bigWig_path, full.names=T, pattern="Preferred")

for (i in 1:length(bigWig_files)) {
#for (i in 1:5) {

  print(paste(i, length(bigWig_files), bigWig_files[i]))
  this_bigWig <- bigWig_files[i]
  thisTarget <- strsplit(this_bigWig, split="/")[[1]]
  thisTarget <- thisTarget[length(thisTarget)]
  thisTarget <- strsplit(thisTarget, split="_")[[1]][1]
  thisTarget <- paste(thisTarget, "_Preferred", sep="")

  thisExpr <- theExprs[grep(thisTarget, theExprs)]
  if(length(thisExpr)!=1) {
    print(i)
  }

  thisName <- strsplit(thisExpr, split="/")[[1]][length(strsplit(thisExpr, split="/")[[1]])]
  thisName <- gsub(".bed.gz", "", thisName)
  thisBed <- thisExpr
  thisBed <- read.delim(thisBed, header=F, sep="\t", stringsAsFactors=F)


  this_cCRE_bound <- paste(outDir_cCREs_Bound, thisName, ".txt", sep="")
  this_Unbound <- paste(outDir_Unbound_cCREs, thisName, ".txt", sep="")

  this_cCRE_bound <- read.table(this_cCRE_bound, header=F, sep="\t", stringsAsFactors=F)
  this_Unbound <- read.table(this_Unbound, header=F, sep="\t", stringsAsFactors=F)

  theCutoff_upper <- as.numeric(quantile(as.numeric(this_cCRE_bound[,2]), 0.9))
  theCutoff_lower <- as.numeric(quantile(as.numeric(this_cCRE_bound[,2]), 0.1))

  theCutoff <- mean(as.numeric(this_cCRE_bound[,2]))

  this_Unbound_high <- this_Unbound[as.numeric(this_Unbound[,2])>=theCutoff_upper,]
  this_Unbound_low <- this_Unbound[as.numeric(this_Unbound[,2])<theCutoff_lower,]


  bed_outfile <- paste(temp_outDir, "peaks.bed", sep="")
  write.table(thisBed[,1:3], bed_outfile, row.names=F, col.names=F, sep="\t", quote=F)

  this_cCRE_bound_locs <- matrix(nrow=nrow(this_cCRE_bound), ncol=3, data=unlist(strsplit(this_cCRE_bound[,1], split=":|-")), byrow=T)
  ccre_bound_outfile <- paste(temp_outDir, "bound_cCREs.bed", sep="")
  write.table(this_cCRE_bound_locs, ccre_bound_outfile, row.names=F, col.names=F, sep="\t", quote=F)

  this_Unbound_low_locs <- matrix(nrow=nrow(this_Unbound_low), ncol=3, data=unlist(strsplit(this_Unbound_low[,1], split=":|-")), byrow=T)
  ccre_Unbound_low_outfile <- paste(temp_outDir, "unbound_cCREs_low.bed", sep="")
  write.table(this_Unbound_low_locs, ccre_Unbound_low_outfile, row.names=F, col.names=F, sep="\t", quote=F)

  this_Unbound_high_locs <- matrix(nrow=nrow(this_Unbound_high), ncol=3, data=unlist(strsplit(this_Unbound_high[,1], split=":|-")), byrow=T)
  ccre_Unbound_high_outfile <- paste(temp_outDir, "unbound_cCREs_high.bed", sep="")
  write.table(this_Unbound_high_locs, ccre_Unbound_high_outfile, row.names=F, col.names=F, sep="\t", quote=F)

  matrix_output <- paste(temp_outDir, "/matrix.mat.gz", sep="")

  print("Making Matrix...")
  theCommand <- paste("computeMatrix reference-point  --referencePoint center -S ", this_bigWig, " -R ", bed_outfile, " ", ccre_bound_outfile, " ", ccre_Unbound_low_outfile, " ", ccre_Unbound_high_outfile, " -a 50 -b 50 --outFileName ", matrix_output, sep="")
  system(theCommand)

  print("Making Profile...")
  figure_saveFile <- paste(figure_outDir, "/", thisTarget, ".pdf", sep="")
  theCommand <- paste("plotProfile -m ", matrix_output, " -out ", figure_saveFile, " --refPointLabel \"Center\" --regionsLabel peaks cCRE_bound cCRE_unbound_l cCRE_unbound_h --legendLocation center-right")
  system(theCommand)

  print("Making Heatmap...")
  figure_saveFile <- paste(figure_outDir, "/", thisTarget, "_heatmap.pdf", sep="")
  theCommand <- paste("plotHeatmap -m ", matrix_output, " -out ", figure_saveFile, " --refPointLabel \"Center\" --regionsLabel peaks cCRE_bound cCRE_unbound_l cCRE_unbound_h --legendLocation center-right")
  system(theCommand)

}
