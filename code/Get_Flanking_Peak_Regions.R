#!/usr/bin/env R
#Get_Flanking_Peak_Regions.R


################################################################################
#This script takes as argument a bed file produced during the running of
#meme_execute.sh, a location to save an additional resulting bed file,
#and a distance by which to extend the edges of peaks for flanking
#sequence. 
################################################################################


options("scipen"=100)

args <- commandArgs(trailingOnly=T)
theBed <- read.table(args[1], header=F, sep="\t", stringsAsFactors=F)
theBed <- as.data.frame(theBed)
theBed[,1] <- as.character(theBed[,1])
theBed[,2] <- as.numeric(as.character(theBed[,2]))
theBed[,3] <- as.numeric(as.character(theBed[,3]))

theBed <- theBed[order(theBed[,1], theBed[,2]),]



starts_1 <- theBed[,2]-(as.numeric(args[3])+1)
ends_1 <- theBed[,2]-1

starts_2 <- theBed[,3]+1
ends_2 <- theBed[,3]+(as.numeric(args[3])+1)


theBed_1 <- theBed
theBed_1[,2] <- starts_1
theBed_1[,3] <- ends_1


theBed_2 <- theBed
theBed_2[,2] <- starts_2
theBed_2[,3] <- ends_2

finalBed <- rbind(theBed_1, theBed_2)
finalBed <- as.data.frame(finalBed)
finalBed[,2] <- as.numeric(as.character(finalBed[,2]))
finalBed[,3] <- as.numeric(as.character(finalBed[,3]))
finalBed[finalBed[,2]<1,2] <- 1


theWidths <- finalBed[,3]-finalBed[,2]
if(length(theWidths[theWidths<as.numeric(args[3])])>0) { finalBed <- finalBed[-which(theWidths<as.numeric(args[3])),] }


finalBed <- finalBed[order(finalBed[,1], finalBed[,2]),]

write.table(finalBed, args[2], row.names=F, col.names=F, sep="\t", quote=F)
