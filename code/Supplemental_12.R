#!/usr/bin/env R
#Supplemental_12.R

################################################################################
#This script is used to produce Supplemental Figure 12
#
#Run under R 4.1.1
#This script takes as arguments:
# outDir: Path where files are to be written
# theMatrix: Table containing information about which TFs had a peak 1kb upstream
#     of the TSS of a given gene, and the TPM of that gene. Provided for download
#     as Binding_Expression_Matrix.txt
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(ggplot2)
library(GenomicRanges)

options("scipen"=100)


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

makeGraph <- function(theY, theX, saveFile, xLab="", yLab="", mainLab="", theCor=NULL, theCor2=NULL) {
  df <- cbind(as.numeric(theY), as.numeric(theX))
  df <- as.data.frame(df, stringsAsFactors=F)
  df$density <- get_density(df[,1], df[,2], n=30)
  ggplot(df, aes(x=theX, y=theY, color = density)) + theme_classic() +
    geom_point() + labs(x=xLab, y=yLab) +
    theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    geom_smooth(method='lm') + xlim(-3,8) + ylim(-3,9)
  ggsave(saveFile)
  #, title=mainLab, subtitle=paste("rho=", round(theCor, digits=4), ", r=", round(theCor2, digits=4), sep="")
}



fitLM <- function(theTable, saveBase) {

  #First, remove any columns in which no TF has bound any promoter.
  theTable <- theTable[,colSums(theTable)>0]

  #Sample the genes to build a model on 70% and test on 30%.
  set.seed(1)
  theSample <- sample(1:nrow(theTable), floor(0.7*nrow(theTable)), replace=F)

  training_Table <- theTable[theSample,]
  test_Table <- theTable[-theSample,]
  theLM <- lm(TPM ~ . , data=training_Table)

  test_data_frame <- as.data.frame(cbind(test_Table[,2:ncol(test_Table)]))
  colnames(test_data_frame) <- c("numBound")
  thePred_lm <- predict(theLM, test_data_frame)
  theCor <- cor.test(test_Table[,1], thePred_lm, method="spearman")
  theCor2 <- cor.test(test_Table[,1], thePred_lm, method="pearson")
  print(theCor2)

  saveFile_1 <- paste(saveBase, ".pdf", sep="")
  mainLab <- strsplit(saveBase, split="/")[[1]]
  mainLab <- mainLab[length(mainLab)]
  xLab <- "Observed log(TPM)"
  yLab <- "Predicted log(TPM)"
  makeGraph(test_Table[,1], thePred_lm, saveFile_1, xLab=xLab, yLab=yLab, mainLab=mainLab, theCor=theCor$estimate, theCor2=theCor2$estimate)

  saveFile_2 <- paste(saveBase, ".model.txt", sep="")
  save(theLM, file=saveFile_2)


}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
theMatrix <- args[2]


#Read in the binding and expression matrix.
binding_expression_matrix <- read.delim(theMatrix, header=T, sep="\t", stringsAsFactors=F)
binding_expression_matrix <- binding_expression_matrix[!is.na(binding_expression_matrix$TPM),]




theMatrix <- cbind(binding_expression_matrix$TPM, rowSums(as.matrix(binding_expression_matrix[,6:ncol(binding_expression_matrix)])))
theMatrix[,1] <- log(theMatrix[,1]+0.1)
theMatrix <- as.data.frame(theMatrix)
colnames(theMatrix) <- c("TPM", "numBound")

saveBase <- paste(outDir, "Supplemental_12", sep="")
fitLM(theMatrix, saveBase)
saveFile_1 <- paste(saveBase, ".pdf", sep="")




#
