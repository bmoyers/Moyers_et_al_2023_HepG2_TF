#!/usr/bin/env R
#Figure_2B_Supplemental_29.R

################################################################################
#This script is used to produce Figure 2B as well as Supplemental Figure 29.
#
#Run under R 4.0.0
#This script takes as arguments:
# outDir: Path where files are to be written
# modelSizeData: Model data produced by binding expression scripts. Provided
#     for download as finding_best_model_size.rds
# bindingExprData: Model data produced by binding expression scripts. Provided
#     for download as binding_expr_models_results.rds
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################


library(ggplot2)


################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
modelSizeData <- args[2]
bindingExprData <- args[3]

theData <- readRDS(modelSizeData)


finalDF <- c()
for (i in 1:length(theData)) {
  thisDF <- theData[[i]]
  thisDF <- cbind(thisDF[,c("pearson_cor", "spearman_cor", "rsq", "rmse", "num_tf")])
  finalDF <- rbind(finalDF, thisDF)
}

finalDF <- as.data.frame(finalDF)
for (i in 1:ncol(finalDF)) { finalDF[,i] <- as.numeric(as.character(finalDF[,i]))}
finalDF$label <- factor(finalDF[,5], levels=seq(5, 200, by=5))

saveFile_1 <- paste(outDir, "Supplemental_29.pdf", sep="")
#png(filename = saveFile_1, width = 1800, height = 1000, units = "px", res = 150, type="cairo")
pdf(saveFile_1)
ggplot(data=finalDF[finalDF[,5]<=100,], aes(x=label, y=pearson_cor)) + theme_classic() +
  geom_boxplot() + labs(x = "Number of predictors", y = "Pearson Correlation") +
    theme(axis.text= element_text(size=20), axis.title=element_text(size=25), axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5))
dev.off()



################
#Okay, now do the actual binding expression results.
################

theData <- readRDS(bindingExprData)

models_df <- as.data.frame(theData[[1]])
summary_df <- as.data.frame(theData[[2]])


graphDF <- c()
for (i in 2:nrow(summary_df)) {
  print(paste(i, nrow(summary_df)))
  thisTF <- summary_df[i,1]
  this_fraction_significant <- summary_df[i,"frac_sig"]
  this_miniModels <- models_df[models_df[,1]==thisTF,c("estimate", "pearson_cor", "spearman_cor", "rsq", "p.value")]
  this_miniModels <- this_miniModels[as.numeric(as.character(this_miniModels[,"p.value"]))<=0.05,]
  median_estimate <- median(as.numeric(as.character(this_miniModels[,"estimate"])))
  theName <- toupper(thisTF)


  miniDF <- cbind(this_miniModels, rep(theName, nrow(this_miniModels)), rep(this_fraction_significant, nrow(this_miniModels)), rep(median_estimate, nrow(this_miniModels)))
  colnames(miniDF) <- c("estimate", "pearson_cor", "spearman_cor", "rsq", "p.value", "TF", "frac_sig", "median_estimate")
  graphDF <- rbind(graphDF, miniDF)
}

graphDF <- as.data.frame(graphDF)
for (i in c(1,2,3,4,5,7,8)) {graphDF[,i] <- as.numeric(as.character(graphDF[,i]))}
graphDF[,6] <- as.character(graphDF[,6])
graphDF[graphDF$TF=="NRSF","TF"] <- "REST"
graphDF[graphDF$TF=="CSDA","TF"] <- "YBX3"


graphDF_50sig <- graphDF[graphDF$frac_sig>=0.5,]

required_top <- c("ATF1", "CEBPG", "FOSL2", "MYC", "NRF1", "ZNF317", "ZNF331")

graphDF_50sig <- graphDF_50sig[order(graphDF_50sig[,"median_estimate"], decreasing=T),]
uniqueTFs <- unique(graphDF_50sig$TF)

which(uniqueTFs%in%required_top)
the_top <- c(uniqueTFs[1:(15-length(required_top))], required_top)

the_top_graphDF <- graphDF_50sig[graphDF_50sig$TF%in%the_top,]



#required_bottom <- c("AHR", "NRSF", "RREB1", "ZFP14", "ZNF660")
required_bottom <- c("AHR", "REST", "RREB1", "ZFP14", "ZNF660")


graphDF_50sig <- graphDF_50sig[order(graphDF_50sig[,"median_estimate"], decreasing=F),]
uniqueTFs <- unique(graphDF_50sig$TF)

which(uniqueTFs%in%required_bottom)
the_bottom <- c(uniqueTFs[1:(15-length(required_bottom))], required_bottom)
the_bottom_graphDF <- graphDF_50sig[graphDF_50sig$TF%in%the_bottom,]




graphDF_combined <- rbind(the_top_graphDF, the_bottom_graphDF)
graphDF_combined <- graphDF_combined[order(graphDF_combined[,"median_estimate"], decreasing=F),]
graphDF_combined$Order <- factor(graphDF_combined$TF, levels=unique(graphDF_combined$TF))



saveFile <- paste(outDir, "Moyers_Figure2B.pdf", sep="")
p <- ggplot(graphDF_combined, aes(x=Order, y=estimate, fill=frac_sig))  + theme_classic() + xlab("TF") + ylab("Model Estimate") +
    geom_boxplot() + coord_flip() + scale_fill_gradient(low="pink", high="lightblue") + ylim(-3,3) +
    theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25,face="bold"), legend.key.size=unit(1, 'cm'), legend.title=element_text(size=25), legend.text=element_text(size=20))
ggsave(saveFile, width=10)






#
