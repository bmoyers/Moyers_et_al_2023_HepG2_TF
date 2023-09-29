#!/usr/bin/env R
#Figure_3A_Supplemental_28.R

################################################################################
#This script is used to produce Figure 3A as well as Supplemental Figure 28.
#
#Run under R 3.6.1
#This script takes as arguments:
# outDir: Path where files are to be written
# abc: Path to ABC connections, provided by Jessie Engreitz
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# finalAnnotationsTFs: Supplemental Table 1
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
library(wesanderson)



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



GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}




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
abc <- args[2]
exprDir <- args[3]
finalAnnotationsTFs <- args[4]

pal <- wes_palette("Moonrise3")


#################################
#First, load in ABC data, and confirm that none of the regions contain
#their own promoters.
#################################


abc <- read.delim(abc, header=T, sep="\t", stringsAsFactors=F)

abc <- abc[as.numeric(as.character(abc$distance))>=10000,]



abc[,1] <- as.character(abc[,1])
abc[,2] <- as.numeric(as.character(abc[,2]))
abc[,3] <- as.numeric(as.character(abc[,3]))

abc_gr <- makeGRangesFromDataFrame(abc, ignore.strand=T, keep.extra.columns=F)





TSS_regions <- cbind(abc[,1], abc[,"TargetGeneTSS"], abc[,"TargetGeneTSS"])
colnames(TSS_regions) <- c("chr", "start", "end")
TSS_regions <- as.data.frame(TSS_regions)
TSS_regions[,2] <- as.numeric(as.character(TSS_regions[,2]))-550
TSS_regions[,3] <- as.numeric(as.character(TSS_regions[,3]))+550
TSS_regions_gr <- makeGRangesFromDataFrame(TSS_regions, ignore.strand=T, keep.extra.columns=T)




#################################
#USe genomic ranges to determine the number of factors bound to each.
#################################

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





#################################
#Look at the overall distribution of factors bound at enhancers-- is it bimodal
#################################


enh_peakNum_perBP <- enh_peakNum/(abc[,3]-abc[,2]+1)
prom_peakNum_perBP <- prom_peakNum/(TSS_regions[,3]-TSS_regions[,2]+1)





################################################################################
#Make a Dataframe of this information, including peak strength and whether or not
#a region is HOT.
################################################################################


theCutoff <- floor(length(peaks)/4)


enh_HOT <- rep("NOT_enh", length(enh_peakNum))
enh_HOT[enh_peakNum>theCutoff] <- "HOT_enh"
prom_HOT <- rep("NOT_prom", length(prom_peakNum))
prom_HOT[prom_peakNum>theCutoff] <- "HOT_prom"

data_frame_abc <- cbind(abc[,1:4], abc[,"TargetGene"], abc[,"ABC.Score"], enh_peakNum, enh_peakNum_perBP, prom_peakNum, prom_peakNum_perBP, prom_perc_inEnh, enh_perc_inProm, jaccard_between, jaccard_between_TF_only, enh_HOT, prom_HOT)
colnames(data_frame_abc) <- c("chr", "start", "end", "name", "gene", "ABC.Score", "enh_peakNum", "enh_peakNum_perBP", "prom_peakNum", "prom_peakNum_perBP", "prom_perc_inEnh", "enh_perc_inProm", "jaccard_between", "jaccard_between_TF_only", "enh_HOT", "prom_HOT")
data_frame_abc <- as.data.frame(data_frame_abc)
data_frame_abc[,1] <- as.character(data_frame_abc[,1])
data_frame_abc[,4] <- as.character(data_frame_abc[,4])
data_frame_abc[,5] <- as.character(data_frame_abc[,5])
data_frame_abc[,15] <- factor(as.character(data_frame_abc[,15]), levels=c("HOT_enh", "NOT_enh"))
data_frame_abc[,16] <- factor(as.character(data_frame_abc[,16]), levels=c("HOT_prom", "NOT_prom"))
for (i in c(2,3,6:14)) { abc[,i] <- as.numeric(as.character(data_frame_abc[,i]))}




HOT <- rep("Neither", nrow(data_frame_abc))
enh_hot <- which(data_frame_abc[,"enh_HOT"]=="HOT_enh")
prom_hot <- which(data_frame_abc[,"prom_HOT"]=="HOT_prom")
both_hot <- enh_hot[enh_hot%in%prom_hot]

#HOT[enh_hot] <- "Partial"
#HOT[prom_hot] <- "Partial"
HOT[enh_hot] <- "Enhancer"
HOT[prom_hot] <- "Promoter"
HOT[both_hot] <- "Both"

data_frame_abc$HOT <- HOT





#################################
#Also look at # of connections per enhancer in HOT enhancers versus NOT.
#################################


data_frame_abc_hotEnh <- data_frame_abc[data_frame_abc[,"enh_HOT"]=="HOT_enh",]
data_frame_abc_notEnh <- data_frame_abc[data_frame_abc[,"enh_HOT"]=="NOT_enh",]
length(unique(data_frame_abc_hotEnh[,"name"]))/nrow(data_frame_abc_hotEnh)
length(unique(data_frame_abc_notEnh[,"name"]))/nrow(data_frame_abc_notEnh)

data_frame_abc_hotProm <- data_frame_abc[data_frame_abc[,"prom_HOT"]=="HOT_prom",]
data_frame_abc_notProm <- data_frame_abc[data_frame_abc[,"prom_HOT"]=="NOT_prom",]
length(unique(data_frame_abc_hotProm[,"name"]))/nrow(data_frame_abc_hotProm)
length(unique(data_frame_abc_notProm[,"name"]))/nrow(data_frame_abc_notProm)


connections_df <- c()
for (i in unique(data_frame_abc_hotEnh[,"name"])) {
  thisSet <- data_frame_abc_hotEnh[data_frame_abc_hotEnh[,"name"]==i,]
  connections_df <- rbind(connections_df, c(i, "HOT Enhancer", thisSet[1,"enh_peakNum"], nrow(thisSet)))
}
for (i in unique(data_frame_abc_notEnh[,"name"])) {
  thisSet <- data_frame_abc_notEnh[data_frame_abc_notEnh[,"name"]==i,]
  connections_df <- rbind(connections_df, c(i, "Non-HOT Enhancer", thisSet[1,"enh_peakNum"], nrow(thisSet)))
}
for (i in unique(data_frame_abc_hotProm[,"gene"])) {
  thisSet <- data_frame_abc_hotProm[data_frame_abc_hotProm[,"gene"]==i,]
  connections_df <- rbind(connections_df, c(i, "HOT Promoter", thisSet[1,"prom_peakNum"], nrow(thisSet)))
}
for (i in unique(data_frame_abc_notProm[,"gene"])) {
  thisSet <- data_frame_abc_notProm[data_frame_abc_notProm[,"gene"]==i,]
  connections_df <- rbind(connections_df, c(i, "Non-HOT Promoter", thisSet[1,"prom_peakNum"], nrow(thisSet)))
}

connections_df <- as.data.frame(connections_df)
connections_df[,2] <- factor(connections_df[,2], levels=c("HOT Enhancer", "Non-HOT Enhancer", "HOT Promoter", "Non-HOT Promoter"))
connections_df[,3] <- as.numeric(as.character(connections_df[,3]))
connections_df[,4] <- as.numeric(as.character(connections_df[,4]))
colnames(connections_df) <- c("name", "type", "numBound", "connections")


my_comparisons <- list( c("HOT Enhancer", "Non-HOT Enhancer"), c("HOT Promoter", "Non-HOT Promoter"))

saveFile_1 <- paste(outDir, "Moyers_Figure3A.pdf", sep="")
p <- ggplot(connections_df[connections_df$connections<=6,], aes(x=type, y=connections)) + theme_classic() +
  theme(axis.text.x=element_text(size=20, angle=90), axis.text.y=element_text(size=20), axis.title.y=element_text(size=25), axis.title.x=element_blank()) + ylim(0,7) +
  geom_boxplot(position=position_dodge(1)) + stat_compare_means(method="t.test", comparisons = my_comparisons, size=6) + ylab("# of ABC Connections")
ggsave(saveFile_1)



##############################################
#I need to remake the last two plots, the Jaccard index plots,
#with a specific comparison to subsampling.  The idea being
#that maybe the sharedness is expected just based upon the large number
#of factors present in the different location.  We need to show
#that it is higher than expected by random chance.
##############################################



renamed_centered_peaks <- c()
for (i in 1:length(peaks)) {
  theName <- strsplit(peaks[i], split="/")[[1]][length(strsplit(peaks[i], split="/")[[1]])]
  theName <- strsplit(theName, split="_")[[1]][1]
  theName <- gsub("-FLAG", "", theName)
  theName <- gsub("-eGFP", "", theName)
  renamed_centered_peaks <- c(renamed_centered_peaks, theName)

}

randomized_Jaccard <- c()
for (i in 1:nrow(data_frame_abc)) {
  if(i%%1000==0) {print(paste(i, nrow(data_frame_abc)))}
  thisProm_num <- data_frame_abc[i,"prom_peakNum"]
  thisEnh_num <- data_frame_abc[i,"enh_peakNum"]
  this_Jaccard <- data_frame_abc[i,"jaccard_between_TF_only"]

  thisProm_sample <- sample(renamed_centered_peaks, thisProm_num, replace=F)
  thisEnh_sample <- sample(renamed_centered_peaks, thisEnh_num, replace=F)

  thisProm_sample <- thisProm_sample[thisProm_sample%in%annotations[,1]]
  thisEnh_sample <- thisEnh_sample[thisEnh_sample%in%annotations[,1]]


  this_rand_Jaccard <- length(thisEnh_sample[thisEnh_sample%in%thisProm_sample])/(length(thisEnh_sample[thisEnh_sample%in%thisProm_sample]) + length(thisEnh_sample[!thisEnh_sample%in%thisProm_sample]) + length(thisProm_sample[!thisProm_sample%in%thisEnh_sample]))
  randomized_Jaccard <- c(randomized_Jaccard, this_rand_Jaccard)


}

data_frame_abc$randomized_Jaccard_TFOnly <- randomized_Jaccard



a <- cbind(data_frame_abc[,c("HOT", "jaccard_between")], rep("observed", nrow(data_frame_abc)))
b <- cbind(data_frame_abc[,c("HOT", "randomized_Jaccard_TFOnly")], rep("randomized", nrow(data_frame_abc)))

colnames(a) <- c("HOT", "Jaccard", "Case")
colnames(b) <- c("HOT", "Jaccard", "Case")
graphDF <- rbind(a, b)
graphDF <- graphDF[!is.na(graphDF$Jaccard),]

graphDF$Label <- paste(graphDF$HOT, graphDF$Case, sep="_")
graphDF$Label <- factor(graphDF$Label, )

graphDF$HOT <- factor(graphDF$HOT, level=c("Both", "Promoter", "Enhancer", "Neither"))


custom_col <- c("green", "grey")

saveFile <- paste(outDir, "Supplemental_28.pdf", sep="")
split_plot<-ggplot(graphDF, aes(x=HOT, y=Jaccard, fill=Case)) + theme_classic() + scale_fill_manual(name="Category", values = custom_col, na.value="grey50") +
  theme(axis.text.x=element_text(size=20, angle=90), axis.text.y=element_text(size=20), axis.title.y=element_text(size=25), axis.title.x=element_text(size=25), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=20), legend.text=element_text(size=15)) +
  geom_split_violin(trim = TRUE, color=NA)+  geom_boxplot(width = 1.0, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0) + ylim(0,.60) + xlab("HOT") + ylab("Promoter:Enhancer Jaccard")
ggsave(saveFile)


cases <- unique(as.character(graphDF$HOT))

pval_table <- c()
for (i in 1:length(cases)) {
  thisSet <- graphDF[graphDF$HOT==cases[i],]
  thisSet_observed <- thisSet[thisSet$Case=="observed",]
  thisSet_random <- thisSet[thisSet$Case=="randomized",]
  this_ttest <- t.test(thisSet_observed$Jaccard, thisSet_random$Jaccard, alternative="two.sided")
  pval_table <- rbind(pval_table, c(cases[i], nrow(thisSet_observed), this_ttest$p.value))
}

colnames(pval_table) <- c("Case", "nObs", "pvalue")
saveFile <- paste(outDir, "Supplemental_Table_17.txt", sep="")
write.table(pval_table, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
