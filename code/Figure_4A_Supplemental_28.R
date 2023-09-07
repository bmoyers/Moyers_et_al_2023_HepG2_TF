#!/usr/bin/env R
#Figure_4A_Supplemental_28.R

################################################################################
#This script is used to produce Figure 4A as well as Supplemental Figure 28.
#
#Run under R 3.6.1
#This script takes as arguments:
# outDir: Path where files are to be written
# noncCRE_table: table of non-cCRE regions with at least one factor bound,
#     provided for download as non_cCRE_allFactors_regions_with_various_annotations.txt
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(ggplot2)
library(matrixStats)
library(ComplexHeatmap)
library(GenomicRanges)

################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
noncCRE_table <- args[2]


#1) Get the table of regions; Restrict to distance >=700bp.
#  Restrict to cases of 2+DAP bound (any DAP excluding histones).

noncCRE_table <- read.delim(noncCRE_table, header=T, sep="\t", stringsAsFactors=F)

noncCRE_table$numFactors <- as.numeric(as.character(noncCRE_table$numFactors))
noncCRE_table$cCRE_dist <- as.numeric(as.character(noncCRE_table$cCRE_dist))

################################################################################
#Make a distribution of distances to cCREs.
################################################################################

saveFile <- paste(outDir, "Supplemental_28.pdf", sep="")
p <- ggplot(noncCRE_table, aes(x=cCRE_dist)) + theme_classic() + geom_density() + xlim(0,10000) + geom_vline(xintercept = 700, color="red")
ggsave(saveFile)


noncCRE_table <- noncCRE_table[noncCRE_table[,"numFactors"]>=2,]
noncCRE_table <- noncCRE_table[noncCRE_table[,"cCRE_dist"]>=700,]


unique(noncCRE_table$other_cCREs)

noncCRE_table <- noncCRE_table[noncCRE_table$other_cCREs%in%c("TF", NA),]
unique(noncCRE_table$other_cCREs)

nrow(noncCRE_table)
#16412




#2) for each region, identify binding of each TF.
#   identify the most common TFs.

theDAPs <- unique(unlist(strsplit(noncCRE_table$Factors, split=",")))

boundTable <- matrix(nrow=nrow(noncCRE_table), ncol=length(theDAPs), data=0)
colnames(boundTable) <- theDAPs

for (i in 1:nrow(noncCRE_table)) {
  thisSet <- unique(unlist(strsplit(noncCRE_table[i,"Factors"], split=",")))
  boundTable[i,thisSet] <- 1
}

theSums <- colSums(boundTable)
boundTable <- boundTable[,order(theSums, decreasing=T)]
theSums <- colSums(boundTable)
length(theSums[theSums>=400])

interesting_boundTable <- boundTable[,theSums>=400]

sMAF_nums <- rowSums(boundTable[,c("MAFF", "MAFK", "MAFG")])


noncCRE_table_multiSMAF <- noncCRE_table[sMAF_nums>=2,]

saveFile <- paste(outDir, "Supplemental_28_multi_sMAF_over700bp.pdf", sep="")
p <- ggplot(noncCRE_table_multiSMAF, aes(x=cCRE_dist)) + theme_classic() + geom_density() + xlim(0,10000) + geom_vline(xintercept = 700, color="red")
ggsave(saveFile)


saveFile <- paste(outDir, "Supplemental_28_multi_sMAF_over700bp_hist.pdf", sep="")
p <- ggplot(noncCRE_table_multiSMAF, aes(x=cCRE_dist)) + theme_classic() + geom_histogram() + xlim(0,10000) + geom_vline(xintercept = 700, color="red")
ggsave(saveFile)

mean(noncCRE_table_multiSMAF$cCRE_dist)



#Create the relevant binding matrix.

heatmapTable <- as.data.frame(interesting_boundTable)
heatmapTable$ATAC <- noncCRE_table$atac_overlap
heatmapTable$H3K9me3 <- noncCRE_table$H3K9me3_overlap
heatmapTable$H3K4me1 <- noncCRE_table$H3K4me1_overlap
heatmapTable$H3K27me3 <- noncCRE_table$H3K27me3_overlap
heatmapTable$H3K27ac <- noncCRE_table$H3K27ac_overlap

RegulatoryMark <- rep(0, nrow(heatmapTable))
RegulatoryMark[heatmapTable$ATAC==1] <- 1
RegulatoryMark[heatmapTable$H3K9me3==1] <- 1
RegulatoryMark[heatmapTable$H3K4me1==1] <- 1
RegulatoryMark[heatmapTable$H3K27me3==1] <- 1
RegulatoryMark[heatmapTable$H3K27ac==1] <- 1

heatmapTable$RegulatoryMark <- RegulatoryMark



heatmapTable <- heatmapTable[,c(colnames(heatmapTable)[1:17], "RegulatoryMark", "H3K9me3", "H3K4me1", "H3K27me3", "H3K27ac", "ATAC")]

heatmapTable <- heatmapTable[order(-heatmapTable$RegulatoryMark, heatmapTable$H3K9me3, heatmapTable$H3K4me1, heatmapTable$ATAC, heatmapTable$H3K27me3, heatmapTable$H3K27ac, heatmapTable$MAFK,
  heatmapTable$MAFF, heatmapTable$ATF7, heatmapTable$JUND, heatmapTable$CREB1, heatmapTable$FOXA2, heatmapTable$HLF, heatmapTable$ONECUT1, heatmapTable$CEBPA, heatmapTable$CEBPB,
  heatmapTable$HNF4A, heatmapTable$ZNF770, heatmapTable$HNRNPK, heatmapTable$ZNF121, heatmapTable$CEBPG, heatmapTable$RBM22, heatmapTable$NFIL3, decreasing=T),]
heatmapTable_matrix <- as.matrix(heatmapTable)

alt_heatmapTable <- heatmapTable
for (i in 1:ncol(alt_heatmapTable)) {
  alt_heatmapTable[heatmapTable[,i]==0,i] <- "Unbound"
  alt_heatmapTable[heatmapTable[,i]==1,i] <- "Bound"

}
heatmapTable_matrix <- as.matrix(alt_heatmapTable)
heatmapTable_matrix <- heatmapTable_matrix[,-grep("RegulatoryMark", colnames(heatmapTable_matrix))]

for (i in 18:21) {
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Unbound",i] <- "No Mark"
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Bound",i] <- "Mark"
}

for (i in 22) {
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Unbound",i] <- "No ATAC"
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Bound",i] <- "ATAC"
}

colors = structure(c("black", "white", "red", "white", "blue", "white"), names = c("Bound", "Unbound", "Mark", "No Mark", "ATAC", "No ATAC"))

this_save_file <- paste(outDir, "Moyers_Figure4A.pdf", sep="")
pdf(this_save_file)
#Heatmap(heatmapTable_matrix, name="Bound",
Heatmap(heatmapTable_matrix, name="Mark in Region",
    column_title="Factor",
    row_title="Regions",
    row_title_gp = gpar(fontsize = 25),
    column_title_gp = gpar(fontsize = 25),
    row_names_gp = gpar(fontsize = 0),
    column_names_side = "top",
    row_names_side= "left",
    column_names_gp = gpar(fontsize = 20),
		cluster_rows=FALSE, cluster_columns=FALSE,
    col = colors, use_raster=FALSE,
    heatmap_legend_param = list(fontsize=15)
)
dev.off()

#Sort things on the side.
#Save using ggplot2 heatmap, I think.


################################################################################
#It has been requested that we also look in to repeat overlaps.
################################################################################


repeats <- read.delim(repeats, header=F, sep="\t", stringsAsFactors=F)
colnames(repeats)[1:3] <- c("chr", "start", "end")
repeats_gr <- makeGRangesFromDataFrame(repeats)
sum(width(reduce(repeats_gr)))
#1636273818
#These repeat regions cover half the genome.

noncCRE_table_gr <- makeGRangesFromDataFrame(noncCRE_table)

repeatOverlap <- rep(0, nrow(noncCRE_table))

theIntersect <- as.data.frame(findOverlaps(noncCRE_table_gr, repeats_gr))
repeatOverlap[unique(theIntersect[,1])] <- 1





heatmapTable <- as.data.frame(interesting_boundTable)
heatmapTable$ATAC <- noncCRE_table$atac_overlap
heatmapTable$H3K9me3 <- noncCRE_table$H3K9me3_overlap
heatmapTable$H3K4me1 <- noncCRE_table$H3K4me1_overlap
heatmapTable$H3K27me3 <- noncCRE_table$H3K27me3_overlap
heatmapTable$H3K27ac <- noncCRE_table$H3K27ac_overlap
heatmapTable$Repeats <- repeatOverlap

RegulatoryMark <- rep(0, nrow(heatmapTable))
RegulatoryMark[heatmapTable$ATAC==1] <- 1
RegulatoryMark[heatmapTable$H3K9me3==1] <- 1
RegulatoryMark[heatmapTable$H3K4me1==1] <- 1
RegulatoryMark[heatmapTable$H3K27me3==1] <- 1
RegulatoryMark[heatmapTable$H3K27ac==1] <- 1
RegulatoryMark[heatmapTable$Repeats==1] <- 1

heatmapTable$RegulatoryMark <- RegulatoryMark



heatmapTable <- heatmapTable[,c(colnames(heatmapTable)[1:17], "RegulatoryMark", "H3K9me3", "H3K4me1", "H3K27me3", "H3K27ac", "ATAC", "Repeats")]

heatmapTable <- heatmapTable[order(heatmapTable$Repeats, -heatmapTable$RegulatoryMark, heatmapTable$H3K9me3, heatmapTable$H3K4me1, heatmapTable$ATAC, heatmapTable$H3K27me3, heatmapTable$H3K27ac, heatmapTable$MAFK,
  heatmapTable$MAFF, heatmapTable$ATF7, heatmapTable$JUND, heatmapTable$CREB1, heatmapTable$FOXA2, heatmapTable$HLF, heatmapTable$ONECUT1, heatmapTable$CEBPA, heatmapTable$CEBPB,
  heatmapTable$HNF4A, heatmapTable$ZNF770, heatmapTable$HNRNPK, heatmapTable$ZNF121, heatmapTable$CEBPG, heatmapTable$RBM22, heatmapTable$NFIL3, decreasing=T),]
heatmapTable_matrix <- as.matrix(heatmapTable)

alt_heatmapTable <- heatmapTable
for (i in 1:ncol(alt_heatmapTable)) {
  alt_heatmapTable[heatmapTable[,i]==0,i] <- "Unbound"
  alt_heatmapTable[heatmapTable[,i]==1,i] <- "Bound"

}
heatmapTable_matrix <- as.matrix(alt_heatmapTable)
heatmapTable_matrix <- heatmapTable_matrix[,-grep("RegulatoryMark", colnames(heatmapTable_matrix))]

for (i in 18:21) {
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Unbound",i] <- "No Mark"
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Bound",i] <- "Mark"
}

for (i in 22) {
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Unbound",i] <- "No ATAC"
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Bound",i] <- "ATAC"
}

for (i in 23) {
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Unbound",i] <- "No Repeat"
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Bound",i] <- "Repeat"
}

colors = structure(c("black", "white", "red", "white", "blue", "white", "green", "white"), names = c("Bound", "Unbound", "Mark", "No Mark", "ATAC", "No ATAC", "Repeat", "No Repeat"))

this_save_file <- paste(outDir, "Moyers_Figure4A_withRepeats.pdf", sep="")
pdf(this_save_file)
#Heatmap(heatmapTable_matrix, name="Bound",
Heatmap(heatmapTable_matrix, name="Mark in Region",
    column_title="Factor",
    row_title="Regions",
    row_title_gp = gpar(fontsize = 25),
    column_title_gp = gpar(fontsize = 25),
    row_names_gp = gpar(fontsize = 0),
    column_names_side = "top",
    row_names_side= "left",
    column_names_gp = gpar(fontsize = 20),
		cluster_rows=FALSE, cluster_columns=FALSE,
    col = colors, use_raster=FALSE,
    heatmap_legend_param = list(fontsize=15)
)
dev.off()



################################################################################
#Compare this to dELS in cCRE.
################################################################################

cCREs <- read.delim(cCREs, header=F, sep="\t", stringsAsFactors=F)

theStates <- unique(cCREs[,6])
repeat_cCRE_overlap_table <- c()
for (i in 1:length(theStates)) {
  thisSet <- cCREs[cCREs[,6]==theStates[i],]
  colnames(thisSet)[1:3] <- c("chr", "start", "end")
  thisSet_gr <- makeGRangesFromDataFrame(thisSet)
  repeatOverlap_thisSet <- rep(0, nrow(thisSet))
  theIntersect <- as.data.frame(findOverlaps(thisSet_gr, repeats_gr))
  repeatOverlap_thisSet[unique(theIntersect[,1])] <- 1
  thisLine <- c(theStates[i], length(repeatOverlap_thisSet), sum(repeatOverlap_thisSet))
  repeat_cCRE_overlap_table <- rbind(repeat_cCRE_overlap_table, thisLine)
}

repeat_cCRE_overlap_table <- as.data.frame(repeat_cCRE_overlap_table)

colnames(repeat_cCRE_overlap_table) <- c("cCRE_State", "numObs", "numOverlapRepeat")
rownames(repeat_cCRE_overlap_table) <- c()
repeat_cCRE_overlap_table$numObs <- as.numeric(as.character(repeat_cCRE_overlap_table$numObs))
repeat_cCRE_overlap_table$numOverlapRepeat <- as.numeric(as.character(repeat_cCRE_overlap_table$numOverlapRepeat))

repeat_cCRE_overlap_table$fractionOverlap <- repeat_cCRE_overlap_table$numOverlapRepeat / repeat_cCRE_overlap_table$numObs







#
