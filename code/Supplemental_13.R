#!/usr/bin/env R
#Supplemental_13.R

################################################################################
#This script is used to produce Supplemental Figure 13.
#
#Run under R 3.6.1
#This script takes as arguments:
# outDir: Path where files are to be written
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# bindingExprData: Model data produced by binding expression scripts. Provided
#     for download as binding_expr_models_results.rds
# cCREs: Path to the Version 4 cCREs provided by Jill Moore in June 2022.
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(GenomicRanges)
library(ggplot2)
library(ggpubr)


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


blank_theme <- theme_minimal()+
  theme(
  axis.text.x=element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold"),
  )


find_TF_name <- function(fileName, theSplit="_") {
  fileName <- gsub(".bed", "", fileName)
  fileName <- gsub("/", "", fileName)
  fileName <- strsplit(fileName, split=theSplit)[[1]]
  #fileName <- paste(fileName[1], fileName[length(fileName)], sep="_")
  return(fileName[1])
}


################################################################################
################################################################################
#Being Script
################################################################################
################################################################################


cCREs <- read.table(cCREs, header=F, sep="\t", stringsAsFactors=F)
cCREs <- cCREs[,c(1:3,6)]
colnames(cCREs)[1:4] <- c("chr", "start", "end", "state")
cCREs_gr <- makeGRangesFromDataFrame(cCREs, ignore.strand=TRUE, keep.extra.columns=F)


theStates <- c("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA-TF", "TF", "CA", "None")
custom_col <- c("#F51313", "#F59913", "#F5C813", "#EEA4A0", "#43B6F3", "#56D59C", "#6EA38B", "#808080", "#000000")


################################################################################
#Read in model and identify the candidate repressors and activators
################################################################################


bindingExprData <- readRDS("/gpfs/gpfs1/home/bmoyers/Figures/ENCODE_500plus/Results_post_July5/BindingExpression/binding_expr_models_results.rds")

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

candidateRepressors <- candidateDF[,6:8]
candidateRepressors <- unique(candidateRepressors)
candidateRepressors <- candidateRepressors[candidateRepressors[,3]<0,]
candidateRepressors <- candidateRepressors[candidateRepressors[,2]>=0.5,]
candidateRepressors <- candidateRepressors[!is.na(candidateRepressors[,1]),]

Repressors <- candidateRepressors[,1]



candidateActivators <- candidateDF[,6:8]
candidateActivators <- unique(candidateActivators)
candidateActivators <- candidateActivators[candidateActivators[,3]>0,]
candidateActivators <- candidateActivators[candidateActivators[,2]>=0.5,]
candidateActivators <- candidateActivators[!is.na(candidateActivators[,1]),]
candidateActivators <- candidateActivators[order(candidateActivators$frac_sig, candidateActivators$median_estimate, decreasing=T),]

Activators <- candidateActivators[,1]
Activators <- Activators[1:26]




################################################################################
#Read in cCREs.
################################################################################


################################################################################
#For all of the Repressors, get the number of peaks in each category.
#Make a circle plot for each.
################################################################################

allPeaks <- list.files(exprDir, full.names=T, pattern="Preferred")

outDir_Indiv <- paste(outDir, "cCRE_distribution_Repressors/", sep="")
if(!file.exists(outDir_Indiv)) {system(paste("mkdir ", outDir_Indiv, sep=""))}

total_peaks_df <- matrix(nrow=length(Repressors), ncol=length(theStates), data=NA)
colnames(total_peaks_df) <- theStates
fraction_peaks_df <- total_peaks_df

for (i in 1:length(Repressors)) {
  thisPeaks <- allPeaks[grep(Repressors[i], allPeaks)]

  this_peaks <- read.table(thisPeaks)
  this_starts <- as.numeric(as.character(this_peaks[,2]))
  this_ends <- as.numeric(as.character(this_peaks[,3]))
  this_means <- (this_starts + this_ends)/2
  this_peaks[,2] <- this_means
  this_peaks[,3] <- this_means
  this_peaks <- as.data.frame(this_peaks)
  colnames(this_peaks)[1:3] <- c("chr", "start", "end")
  this_peaks_gr <- makeGRangesFromDataFrame(this_peaks, ignore.strand=T)

  cCREs_Intersect <- as.data.frame(findOverlaps(this_peaks_gr, cCREs_gr))
  these_states <- cCREs[cCREs_Intersect[,2],4]
  these_states_counts <- c()
  for (j in 1:length(theStates)) {
    these_states_counts[j] <- length(these_states[these_states==theStates[j]])
  }
  these_states_counts[length(these_states_counts)] <- nrow(this_peaks)-sum(these_states_counts)

  total_peaks_df[i,] <- these_states_counts
  these_states_counts_fractions <- these_states_counts/nrow(this_peaks)
  fraction_peaks_df[i,] <- these_states_counts_fractions

  this_save_file <- paste(outDir_Indiv, Repressors[i], "_cCREs_piechart.pdf", sep="")
  piechart_label <- paste(Repressors[i], " distribution with cCREs labels", sep="")
  piechart_sub <- paste("In ", nrow(this_peaks), " peaks", sep="")
  df <- data.frame(cCREs=theStates, Fraction=as.numeric(these_states_counts_fractions))
	df$cCREs <- factor(df$cCREs, levels=df$cCREs)
	bp <- ggplot(df, aes(x="", y=Fraction, fill=cCREs)) + geom_bar(stat="identity", color="black", position=position_fill(reverse = TRUE)) + scale_fill_manual(name="cCREs", values = custom_col, na.value="grey50")
	pie <- bp + coord_polar("y", start=0) + labs(title=piechart_label, subtitle=piechart_sub) + blank_theme
  ggsave(this_save_file)


}

rownames(total_peaks_df) <- Repressors
rownames(fraction_peaks_df) <- Repressors

################################################################################
#Make 2 barplots for repressors, one with total peak counts and one with fractions,
#showing cCRE distribution for each factor.
################################################################################


graphDF <- c()
for (i in 1:nrow(total_peaks_df)) {
  for (j in 1:ncol(total_peaks_df)) {
    thisLine <- c(rownames(total_peaks_df)[i], colnames(total_peaks_df)[j], total_peaks_df[i,j], fraction_peaks_df[i,j])
    graphDF <- rbind(graphDF, thisLine)
  }
}
graphDF <- as.data.frame(graphDF)
colnames(graphDF) <- c("Repressor", "cCRE", "PeakCount", "PeakFraction")
graphDF[,3] <- as.numeric(as.character(graphDF[,3]))
graphDF[,4] <- as.numeric(as.character(graphDF[,4]))
graphDF$cCRE <- factor(graphDF$cCRE, levels=theStates)
graphDF$Repressor <- factor(graphDF$Repressor, levels=Repressors)

saveFile <- paste(outDir, "Supplemental_13_Repressors.pdf", sep="")
ggplot(data=graphDF, aes(x=Repressor, y=PeakCount, fill=cCRE)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25), legend.title=element_text(size=15), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Repressor") + ylab("Number of Peaks") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)


saveFile <- paste(outDir, "cCRE_distribution_peakFraction_Repressors.pdf", sep="")
ggplot(data=graphDF, aes(x=Repressor, y=PeakFraction, fill=cCRE)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25), legend.title=element_text(size=15), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Repressor") + ylab("Fraction of Peaks") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)


graphDF_repressors <- graphDF




################################################################################
#For all of the Activators, get the number of peaks in each category.
#Make a circle plot for each.
################################################################################

allPeaks <- list.files(exprDir, full.names=T, pattern="Preferred")

outDir_Indiv <- paste(outDir, "cCRE_distribution_Activators/", sep="")
if(!file.exists(outDir_Indiv)) {system(paste("mkdir ", outDir_Indiv, sep=""))}

total_peaks_df <- matrix(nrow=length(Activators), ncol=length(theStates), data=NA)
colnames(total_peaks_df) <- theStates
fraction_peaks_df <- total_peaks_df

for (i in 1:length(Activators)) {
  thisPeaks <- allPeaks[grep(Activators[i], allPeaks)]
  if(length(thisPeaks)>1) {
    thisPeaks <- thisPeaks[grep(paste("/", Activators[i], sep=""), thisPeaks)]
  }

  this_peaks <- read.table(thisPeaks)
  this_starts <- as.numeric(as.character(this_peaks[,2]))
  this_ends <- as.numeric(as.character(this_peaks[,3]))
  this_means <- (this_starts + this_ends)/2
  this_peaks[,2] <- this_means
  this_peaks[,3] <- this_means
  this_peaks <- as.data.frame(this_peaks)
  colnames(this_peaks)[1:3] <- c("chr", "start", "end")
  this_peaks_gr <- makeGRangesFromDataFrame(this_peaks, ignore.strand=T)

  cCREs_Intersect <- as.data.frame(findOverlaps(this_peaks_gr, cCREs_gr))
  these_states <- cCREs[cCREs_Intersect[,2],4]
  these_states_counts <- c()
  for (j in 1:length(theStates)) {
    these_states_counts[j] <- length(these_states[these_states==theStates[j]])
  }
  these_states_counts[length(these_states_counts)] <- nrow(this_peaks)-sum(these_states_counts)

  total_peaks_df[i,] <- these_states_counts
  these_states_counts_fractions <- these_states_counts/nrow(this_peaks)
  fraction_peaks_df[i,] <- these_states_counts_fractions

  this_save_file <- paste(outDir_Indiv, Activators[i], "_cCREs_piechart.pdf", sep="")
  piechart_label <- paste(Activators[i], " distribution with cCREs labels", sep="")
  piechart_sub <- paste("In ", nrow(this_peaks), " peaks", sep="")
  df <- data.frame(cCREs=theStates, Fraction=as.numeric(these_states_counts_fractions))
	df$cCREs <- factor(df$cCREs, levels=df$cCREs)
	bp <- ggplot(df, aes(x="", y=Fraction, fill=cCREs)) + geom_bar(stat="identity", color="black", position=position_fill(reverse = TRUE)) + scale_fill_manual(name="cCREs", values = custom_col, na.value="grey50")
	pie <- bp + coord_polar("y", start=0) + labs(title=piechart_label, subtitle=piechart_sub) + blank_theme
  ggsave(this_save_file)


}

rownames(total_peaks_df) <- Activators
rownames(fraction_peaks_df) <- Activators


############################################################################
#Make 2 barplots for activators, one with total peak counts and one with fractions,
#showing cCRE distribution for each factor.
################################################################################


graphDF <- c()
for (i in 1:nrow(total_peaks_df)) {
  for (j in 1:ncol(total_peaks_df)) {
    thisLine <- c(rownames(total_peaks_df)[i], colnames(total_peaks_df)[j], total_peaks_df[i,j], fraction_peaks_df[i,j])
    graphDF <- rbind(graphDF, thisLine)
  }
}
graphDF <- as.data.frame(graphDF)
colnames(graphDF) <- c("Activator", "cCRE", "PeakCount", "PeakFraction")
graphDF[,3] <- as.numeric(as.character(graphDF[,3]))
graphDF[,4] <- as.numeric(as.character(graphDF[,4]))
graphDF$cCRE <- factor(graphDF$cCRE, levels=theStates)
graphDF$Activator <- factor(graphDF$Activator, levels=Activators)

saveFile <- paste(outDir, "Supplemental_13_Activators.pdf", sep="")
ggplot(data=graphDF, aes(x=Activator, y=PeakCount, fill=cCRE)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25), legend.title=element_text(size=15), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Activator") + ylab("Number of Peaks") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)


saveFile <- paste(outDir, "cCRE_distribution_peakFraction_Activators.pdf", sep="")
ggplot(data=graphDF, aes(x=Activator, y=PeakFraction, fill=cCRE)) + geom_bar(stat="identity") + theme_classic() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=25), legend.title=element_text(size=15), legend.text=element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Activator") + ylab("Number of Peaks") + scale_fill_manual(name="cCRE class", values = custom_col, na.value="grey50")
ggsave(saveFile)


graphDF_activators <- graphDF


############################################################################
#Do a fraction comparison of all Repressors and Activators in each cCRE type.
############################################################################


graphDF_cCRE_fractions <- c()

for (i in 1:length(theStates)) {
  these_repressors <- graphDF_repressors[graphDF_repressors$cCRE==theStates[i],c("PeakFraction", "cCRE")]
  these_repressors <- as.data.frame(cbind(these_repressors, rep("Repressor", nrow(these_repressors))))
  colnames(these_repressors)[3] <- "Type"
  graphDF_cCRE_fractions <- rbind(graphDF_cCRE_fractions, these_repressors)

  these_activators <- graphDF_activators[graphDF_activators$cCRE==theStates[i],c("PeakFraction", "cCRE")]
  these_activators <- as.data.frame(cbind(these_activators, rep("Activator", nrow(these_activators))))
  colnames(these_activators)[3] <- "Type"
  graphDF_cCRE_fractions <- rbind(graphDF_cCRE_fractions, these_activators)

}



saveFile <- paste(outDir, "cCRE_distribution_peakFraction_TypeComparison_boxplot.pdf", sep="")
p <- ggplot(graphDF_cCRE_fractions, aes(x=cCRE, y=PeakFraction, fill=Type)) + theme_classic() + geom_boxplot() + ylab("PeakFraction") + xlab("cCREs") +
  theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.title=element_text(size=25), legend.text=element_text(size=20)) + ylim(0,0.75) +
  stat_compare_means(aes(group=Type), label = "p.signif", label.y = 0.7, method="wilcox.test")
ggsave(saveFile)
