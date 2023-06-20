#Fig4_A.R


################################################################################
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



heatmapTable <- heatmapTable[,c(colnames(heatmapTable)[1:17], "RegulatoryMark", "H3K9me3", "H3K4me1", "ATAC", "H3K27me3", "H3K27ac")]

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

for (i in 18:22) {
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Unbound",i] <- "No Mark"
  heatmapTable_matrix[heatmapTable_matrix[,i]=="Bound",i] <- "Mark"

}

colors = structure(c("darkgreen", "white", "Yellow", "white"), names = c("Bound", "Unbound", "Mark", "No Mark")) # Blue, Red

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


#
