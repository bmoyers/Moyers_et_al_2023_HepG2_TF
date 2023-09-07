#!/usr/bin/env R
#Supplemental_11.R

################################################################################
#This script is used to produce Supplemental Figure 11
#
#Run under R 4.1.1
#This script takes as arguments:
# outDir: Path where files are to be written
# thePromoters: file containing refseq TSS with gene names +/- 500bp,
#     provided for download as refseq_genes_unique_TSS_1000.bed
# expr1 and expr2: Path to expression levels of genes in HepG2,
#     accession numbers ENCFF533XPJ and ENCFF321JIT
# expr1_k562 and expr2_k562: Path to expression levels of genes in HepG2,
#     accession numbers ENCFF286GLL and ENCFF986DBN
# finalAnnotationsTFs: Supplemental Table 1
# exprDir: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set.
# exprDir_k562: Path to the directory containing all ChIP-seq experiments in HepG2,
#     provided for download in folder Experiment_Set_K562.
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(GenomicRanges)
library(biomaRt)
library(ggplot2)
library(MASS)
library(matrixStats)

################################################################################
################################################################################
#Functions
################################################################################
################################################################################


################################################################################
#This function takes a dataframe of promoters associated with genes (must have
#gene names in the "name" column), a table of expression values, and a table of genes
#relating gene names to gene expression values.  It then goes through and determines
#the maximum expression level associated with a given promoter.
################################################################################
getNamesAndExpr <- function(thisCounts, thisExpr, theGenes) {
  thisCounts <- read.table(thisCounts, header=F, sep="\t", stringsAsFactors=F)
  colnames(thisCounts) <- c("chr", "start", "end", "name", "IDK", "zstrand")
  expressionVector <- rep(NA, nrow(thisCounts))
  for (i in 1:nrow(thisCounts)) {
    if(i%%1000==0) { print(paste(i, nrow(thisCounts), sep="  "))}
    this_ENSG <- unique(theGenes[which(theGenes[,3]==thisCounts[i,"name"]),1])
    if(length(this_ENSG)>0) {
      exprSet <- c()
      for (j in 1:length(this_ENSG)) {
        exprSet <- rbind(exprSet, thisExpr[grep(this_ENSG[j], thisExpr[,1], fixed=T),])
      }
      if(length(unique(exprSet[,"TPM"]))==1) {
        expressionVector[i] <- exprSet[1,"TPM"]
      }
      if(length(unique(exprSet[,"TPM"]))>1) {
        expressionVector[i] <- max(exprSet[,"TPM"])
      }
    }
    if(length(this_ENSG)==0) {
      expressionVector[i] <- NA
    }
  }
  TPM <- expressionVector
  countsTable <- cbind(thisCounts, TPM)
  return(countsTable)
}



################################################################################
#The following functions make graphs of predicted versus observed data.
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
  ggplot(df, aes(x=theX, y=theY, color = density)) +
    geom_point() + labs(x=xLab, y=yLab) + theme_classic() +
    theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    geom_smooth(method='lm') #+ xlim(-3,8)
  ggsave(saveFile)
  #, title=mainLab, subtitle=paste("rho=", round(theCor, digits=4), ", r=", round(theCor2, digits=4), sep="")
}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
thePromoters <- args[2]
expr1 <- args[3]
expr2 <- args[4]
expr1_k562 <- args[5]
expr2_k562 <- args[6]
finalAnnotationsTFs <- args[7]
exprDir <- args[8]
exprDir_k562 <- args[9]



################################################################################
#Get the TSS+/-1kb file;
#Get expression level for each gene.
################################################################################

TSS <- thePromoters

expr1 <- read.delim(expr1, header=T, sep="\t", stringsAsFactors=F)
expr2 <- read.delim(expr2, header=T, sep="\t", stringsAsFactors=F)
finalExpr <- cbind(expr1[,1:3], c((expr1[,4]+expr2[,4])/2), c((expr1[,5]+expr2[,5])/2), c((expr1[,6]+expr2[,6])/2))
colnames(finalExpr) <- c(colnames(expr1)[1:6])
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
the_genes <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)
TSS_expr <- getNamesAndExpr(TSS, finalExpr, the_genes)


TSS_expr <- TSS_expr[,c(1:4,7)]
TSS_expr <- as.data.frame(TSS_expr)
TSS_expr_gr <- makeGRangesFromDataFrame(TSS_expr, ignore.strand=T, keep.extra.columns=F)

################################################################################
#Identify the relevant files, preferred only.
#remove histones.
#restrict to TFs only.
################################################################################

################################################################################
#HepG2
################################################################################

annotations <- read.delim(finalAnnotationsTFs, header=T, sep="\t", stringsAsFactors=F)
annotations <- annotations[,c(1,3,5)]
annotations <- annotations[annotations[,2]=="Preferred",]
annotations <- annotations[annotations[,3]=="TF",]

experimentFiles <- list.files(exprDir, pattern="Preferred", full.names=T)

theColnames <- c()
theBound <- c()
for (i in 1:length(experimentFiles)) {
  these_bound <- rep(0, nrow(TSS_expr))
  thisTF <- strsplit(experimentFiles[i], split="/")[[1]]
  thisTF <- thisTF[length(thisTF)]
  thisTF <- strsplit(thisTF, split="_")[[1]][1]
  thisTF <- gsub("-FLAG", "", thisTF)
  thisTF <- gsub("-eGFP", "", thisTF)
  if(!thisTF%in%annotations[,1]) {print(thisTF)}
  if(thisTF%in%annotations[,1]) {
    theColnames <- c(theColnames, thisTF)
    thisFile <- read.delim(experimentFiles[i], header=F, sep="\t", stringsAsFactors=F)
    thisFile <- as.data.frame(thisFile)
    colnames(thisFile)[1:3] <- c("chr", "start", "end")
    thisFile_gr <- makeGRangesFromDataFrame(thisFile, ignore.strand=T, keep.extra.columns=F)
    theIntersect <- as.data.frame(findOverlaps(TSS_expr_gr, thisFile_gr))

    these_bound[unique(theIntersect[,1])] <- 1

    theBound <- cbind(theBound, these_bound)
  }

}

colnames(theBound) <- theColnames


final_table <- cbind(TSS_expr, theBound)




################################################################################
#Now do K562
################################################################################

expr1 <- read.delim(expr1_k562, header=T, sep="\t", stringsAsFactors=F)
expr2 <- read.delim(expr1_k562, header=T, sep="\t", stringsAsFactors=F)
finalExpr <- cbind(expr1[,1:3], c((expr1[,4]+expr2[,4])/2), c((expr1[,5]+expr2[,5])/2), c((expr1[,6]+expr2[,6])/2))
colnames(finalExpr) <- c(colnames(expr1)[1:6])
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
the_genes <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)
TSS_expr <- getNamesAndExpr(TSS, finalExpr, the_genes)


TSS_expr <- TSS_expr[,c(1:4,7)]
TSS_expr <- as.data.frame(TSS_expr)
TSS_expr_gr <- makeGRangesFromDataFrame(TSS_expr, ignore.strand=T, keep.extra.columns=F)


experimentFiles <- list.files(exprDir_k562, full.names=T)

theColnames <- c()
theBound <- c()
for (i in 1:length(experimentFiles)) {
  these_bound <- rep(0, nrow(TSS_expr))
  thisTF <- strsplit(experimentFiles[i], split="/")[[1]]
  thisTF <- thisTF[length(thisTF)]
  thisTF <- strsplit(thisTF, split="_")[[1]][1]
  thisTF <- gsub("-FLAG", "", thisTF)
  thisTF <- gsub("-eGFP", "", thisTF)
  if(!thisTF%in%annotations[,1]) {print(thisTF)}
  if(thisTF%in%annotations[,1]) {
    theColnames <- c(theColnames, thisTF)
    thisFile <- read.delim(experimentFiles[i], header=F, sep="\t", stringsAsFactors=F)
    thisFile <- as.data.frame(thisFile)
    colnames(thisFile)[1:3] <- c("chr", "start", "end")
    thisFile_gr <- makeGRangesFromDataFrame(thisFile, ignore.strand=T, keep.extra.columns=F)
    theIntersect <- as.data.frame(findOverlaps(TSS_expr_gr, thisFile_gr))

    these_bound[unique(theIntersect[,1])] <- 1

    theBound <- cbind(theBound, these_bound)
  }

}

colnames(theBound) <- theColnames


final_table_k562 <- cbind(TSS_expr, theBound)


################################################################################
#Restrict each matrix to the same factors.
#train in HepG2
#test in K562.
#To avoid some problem of overfitting, do the same, but use only 70%/30% of the data,
#nonoverlapping, in each cell type.
################################################################################


final_table_hepG2 <- final_table[,colnames(final_table)%in%colnames(final_table_k562)]

final_table_hepG2 <- final_table_hepG2[,5:ncol(final_table_hepG2)]
final_table_k562 <- final_table_k562[,5:ncol(final_table_k562)]

final_table_hepG2$TPM <- log(final_table_hepG2$TPM+0.01)
final_table_hepG2 <- final_table_hepG2[!is.na(final_table_hepG2$TPM),]
final_table_k562$TPM <- log(final_table_k562$TPM+0.01)
final_table_k562 <- final_table_k562[!is.na(final_table_k562$TPM),]


#Using a 70/30 split in the data.
set.seed(1)
theSample <- sample(1:nrow(final_table_hepG2), floor(0.7*nrow(final_table_hepG2)), replace=F)

final_table_hepG2_test <- final_table_hepG2[theSample,]
final_table_k562_pred <- final_table_k562[-theSample,]


theLM <- lm(TPM ~ . , data=final_table_hepG2_test)

thePred_lm <- predict(theLM, final_table_k562_pred[,2:ncol(final_table_k562_pred)])
theCor <- cor.test(final_table_k562_pred[,1], thePred_lm, method="spearman")
theCor2 <- cor.test(final_table_k562_pred[,1], thePred_lm, method="pearson")
print(theCor2)

xLab <- "Observed log(TPM)"
yLab <- "Predicted log(TPM)"

saveFile_1 <- paste(outDir, "Supplemental_11.pdf", sep="")
makeGraph(final_table_k562_pred[,1], thePred_lm, saveFile_1, xLab=xLab, yLab=yLab)


#
#
