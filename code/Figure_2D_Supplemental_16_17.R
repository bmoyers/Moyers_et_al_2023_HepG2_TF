#!/usr/bin/env R
#Figure_2D_Supplemental_16_17.R


################################################################################
#This script is used to produce Figure 2D as well as Supplemental Figures 16 and 17.
#
#Run under R 4.1.1
#This script takes as arguments:
# outDir: Path where files are to be written
# dataDir: Path to directory containing the element representation files,
#     provided for download. Files should begin with the string "Element_representation"
# the_info_table: Supplemental Table 5
# bindingExprData: Model data produced by binding expression scripts. Provided
#     for download as binding_expr_models_results.rds
#
################################################################################


################################################################################
################################################################################
#Load Libraries.
################################################################################
################################################################################

library(rtracklayer)
library(universalmotif)
library(memes)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(plyranges)
library(ggplot2)



################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
dataDir <- args[2]
the_info_table <- args[3]


################################################################################
#Read in the info Table.
################################################################################


the_info_table <- read.table(the_info_table, header=T, sep="\t", stringsAsFactors=F)
the_info_table[,1] <- gsub(">", "", the_info_table[,1])
the_info_table[,2] <- gsub(">", "", the_info_table[,2])

the_info_table_2 <- the_info_table

the_info_table[,"shortTitles"] <- paste(the_info_table[,"shortTitles"], "_fwd", sep="")
the_info_table_2[,"shortTitles"] <- paste(the_info_table_2[,"shortTitles"], "_rev", sep="")
the_info_table_2[,"sequence"] <- as.character(reverseComplement(DNAStringSet(the_info_table_2[,"sequence"])))

the_info_table <- rbind(the_info_table, the_info_table_2)






################################################################################
#We will need to read in predictions based on our binding-expression models
#to identify the candidate repressors and activators.
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

candidateDF <- unique(candidateDF[,6:8])

Repressors <- candidateDF[candidateDF$median_estimate<0,"TF"]
Activators <- candidateDF[candidateDF$median_estimate>0,"TF"]


################################################################################
#Read in the data for individual reps.
################################################################################

################################################################################
#Read in the data for individual reps.
################################################################################


element_saveFiles <- list.files(dataDir, full.names=T, pattern="Element_representation")


element_rep_list <- list()
for (i in 1:length(element_saveFiles)) {
  thisElement <- read.table(element_saveFiles[i], header=T, sep="\t", stringsAsFactors=F)
  element_rep_list[[length(element_rep_list)+1]] <- thisElement
}




################################################################################
#Set up controls
################################################################################

controls <- the_info_table[grep("ontrol", the_info_table[,2]),]
controls_positive <- controls[grep("positive", controls[,2]),1]
controls_negative <- controls[grep("negative", controls[,2]),1]

controls_positive <- gsub(">", "", controls_positive)
controls_positive <- c(paste(controls_positive, "_fwd", sep=""), paste(controls_positive, "_rev", sep=""))

controls_negative <- gsub(">", "", controls_negative)
controls_negative <- c(paste(controls_negative, "_fwd", sep=""), paste(controls_negative, "_rev", sep=""))


combinedTable <- cbind(element_rep_list[[1]][,c("Element", "log2")], element_rep_list[[2]][,"log2"], element_rep_list[[3]][,"log2"])
colnames(combinedTable) <- c("name", "Rep1_log2Ratio", "Rep2_log2Ratio", "Rep3_log2Ratio")
combinedTable <- as.data.frame(combinedTable)
labels <- rep("Element", nrow(combinedTable))
labels[which(combinedTable$name%in%controls_positive)] <- "positive"
labels[which(combinedTable$name%in%controls_negative)] <- "negative"
combinedTable$labels <- labels



################################################################################
#Add data for summing across reps.
################################################################################

dna_counts <- rowSums(cbind(element_rep_list[[1]][,"dna_counts"], element_rep_list[[2]][,"dna_counts"], element_rep_list[[3]][,"dna_counts"]))
rna_counts <- rowSums(cbind(element_rep_list[[1]][,"rna_counts"], element_rep_list[[2]][,"rna_counts"], element_rep_list[[3]][,"rna_counts"]))

combinedMat <- cbind(element_rep_list[[1]][,"Element"], rna_counts, dna_counts)
colnames(combinedMat) <- c("Element", "rna_counts", "dna_counts")
combinedMat <- as.data.frame(combinedMat)
combinedMat$rna_counts <- as.numeric(combinedMat$rna_counts)
combinedMat$dna_counts <- as.numeric(combinedMat$dna_counts)
combinedMat$norm_rna_counts <- combinedMat$rna_counts / (sum(combinedMat$rna_counts)/1000000)
combinedMat$norm_dna_counts <- combinedMat$dna_counts / (sum(combinedMat$dna_counts)/1000000)
combinedMat$ratio <- (combinedMat$norm_rna_counts+0.01) / (combinedMat$norm_dna_counts+0.01)
combinedMat$log2 <- log2(combinedMat$ratio)



combinedTable$sum_log2 <- combinedMat$log2




################################################################################
#Add the category and fullTitles.
################################################################################


category <- c()
for (i in 1:nrow(the_info_table)) {
  thisOne <- strsplit(the_info_table[i,"fullTitles"], split="_")[[1]][2]
  if(thisOne%in%c("H", "ENH")) {
    thisOne <- "bindingExpr"
  }
  if(thisOne=="control") {
    if(length(grep("positive", the_info_table[i,"fullTitles"]))==1) { thisOne <- "positive" }
    if(length(grep("negative", the_info_table[i,"fullTitles"]))==1) { thisOne <- "negative" }
    if(length(grep("bindingExpr", the_info_table[i,"fullTitles"]))==1) { thisOne <- "bindingExpr" }

  }
  category <- c(category, thisOne)
}

category[grep("motifMut", the_info_table$fullTitles)] <- "motifMut"
the_info_table$category <- category


combinedTable$category <- the_info_table[match(combinedTable$name, the_info_table$shortTitles),"category"]
combinedTable$fullTitles <- the_info_table[match(combinedTable$name, the_info_table$shortTitles),"fullTitles"]


################################################################################
#Restrict to the relevant BindingExpression set.
################################################################################


combinedTable_bindingExpr <- combinedTable[grep("bindingExpr", combinedTable$fullTitles),]
combinedTable_bindingExpr_control <- combinedTable_bindingExpr[grep("control", combinedTable_bindingExpr$fullTitles),]
combinedTable_bindingExpr <- combinedTable_bindingExpr[-grep("control", combinedTable_bindingExpr$fullTitles),]





################################################################################
################################################################################
#Simple number analysis.
################################################################################
################################################################################

theTF <- c()
theNumber <- c()
for (i in 1:nrow(combinedTable_bindingExpr)) {
  thisSet <- strsplit(combinedTable_bindingExpr[i,"fullTitles"], split=":", fixed=T)[[1]]
  thisSet <- thisSet[3:length(thisSet)]
  thisSet <- thisSet[seq(1, length(thisSet), by=2)]
  theNumber <- c(theNumber, length(thisSet))
  theTF <- c(theTF, strsplit(thisSet[1], split="_")[[1]][1])

}
uniqueTFs <- unique(theTF)

Repressors <- c(Repressors, "NRSF")
candidateActivators <- uniqueTFs[uniqueTFs%in%Activators]
candidateRepressors <- uniqueTFs[uniqueTFs%in%Repressors]


uniqueTFs <- c(candidateActivators, candidateRepressors)


combinedTable_bindingExpr$TF <- theTF
combinedTable_bindingExpr$numMotifs <- theNumber
combinedTable_bindingExpr_control$TF <- rep("Control", nrow(combinedTable_bindingExpr_control))
combinedTable_bindingExpr_control$numMotifs <- rep(0, nrow(combinedTable_bindingExpr_control))

#
#Make the above a ratio of observation/mean(controls)
#

theMean <- mean(combinedTable_bindingExpr_control$sum_log2)


graphDF <- rbind(combinedTable_bindingExpr)

graphDF$numMotifs <- factor(graphDF$numMotifs, levels=unique(graphDF$numMotifs))
graphDF$TF <- factor(graphDF$TF, levels=c("Control", uniqueTFs))
graphDF$log2_overControl <- graphDF$sum_log2 / theMean



newOrder <- c("FOSL2", "NFYC", "NRF1", "CEBPG", "ATF1", "MYC", "ZNF331", "ZNF317", "NRSF", "ZNF660", "RREB1", "ZFP14", "AHR")
namesColors <- c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "red", "red", "red", "red", "red")
graphDF$TF <- as.character(graphDF$TF)
graphDF$TF <- factor(graphDF$TF, levels=newOrder)

saveFile <- paste(outDir, "Moyers_Figure2D.pdf", sep="")
p <- ggplot(graphDF, aes(x=TF, y=log2_overControl, fill=numMotifs)) + theme_bw() + geom_boxplot() + ylab("Element / Control") + xlab("Motif") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20)) + ylim(-3,6.5) +
  geom_hline(yintercept=1) + scale_fill_grey(start = 0.7, end = 0.3)
ggsave(saveFile)



################################################################################
#Break this down by the particular promoter type, too.
################################################################################

graphDF_ENH <- graphDF[grep("1_binding", graphDF$fullTitles),]

thisMean <- mean(combinedTable_bindingExpr_control[grep("1_binding", combinedTable_bindingExpr_control$fullTitles),"sum_log2"])

graphDF_ENH$log2_overControl_alt <- graphDF_ENH$sum_log2 / thisMean


#
#Make the above a ratio of observation/mean(controls)
#

newOrder <- c("FOSL2", "NFYC", "NRF1", "CEBPG", "ATF1", "MYC", "ZNF331", "ZNF317", "NRSF", "ZNF660", "RREB1", "ZFP14", "AHR")
graphDF_ENH$TF <- as.character(graphDF_ENH$TF)
graphDF_ENH$TF <- factor(graphDF_ENH$TF, levels=newOrder)

saveFile <- paste(outDir, "Supplemental_17.pdf", sep="")
p <- ggplot(graphDF_ENH, aes(x=TF, y=log2_overControl_alt, fill=numMotifs)) + theme_bw() + geom_boxplot() + ylab("MPRA Signal") + xlab("Motif") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, color=namesColors), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20)) + #ylim(-3,6.5) +
  geom_hline(yintercept=1) + scale_fill_grey(start = 0.7, end = 0.3)
ggsave(saveFile)




graphDF_H_046 <- graphDF[grep("046_binding", graphDF$fullTitles),]

thisMean <- mean(combinedTable_bindingExpr_control[grep("046_binding", combinedTable_bindingExpr_control$fullTitles),"sum_log2"])

graphDF_H_046$log2_overControl_alt <- graphDF_H_046$sum_log2 / thisMean

#
#Make the above a ratio of observation/mean(controls)
#




newOrder <- c("FOSL2", "NFYC", "NRF1", "CEBPG", "ATF1", "MYC", "ZNF331", "ZNF317", "NRSF", "ZNF660", "RREB1", "ZFP14", "AHR")
graphDF_H_046$TF <- as.character(graphDF_H_046$TF)
graphDF_H_046$TF <- factor(graphDF_H_046$TF, levels=newOrder)


saveFile <- paste(outDir, "Supplemental_16.pdf", sep="")
p <- ggplot(graphDF_H_046, aes(x=TF, y=log2_overControl_alt, fill=numMotifs)) + theme_bw() + geom_boxplot() +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, color=namesColors), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20)) + #ylim(-3,6.5) +
  geom_hline(yintercept=1) + scale_fill_grey(start = 0.7, end = 0.3)
ggsave(saveFile)





#
