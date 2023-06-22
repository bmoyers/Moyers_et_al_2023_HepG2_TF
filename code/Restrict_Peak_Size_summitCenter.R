#!/usr/bin/env R
#Restrict_Peak_Size_summitCenter.R

################################################################################
#This script takes as argument a bed file produced during the running of
#lsgkm_execute.sh and meme_execute.sh.  It identifies the center of a peak and
#extends the peak by a given number of basepairs.
################################################################################


options("scipen"=100)

args <- commandArgs(trailingOnly=T)
theBed <- read.table(args[1], header=F, sep="\t", stringsAsFactors=F)
theBed <- as.data.frame(theBed)
theBed[,1] <- as.character(theBed[,1])
theBed[,2] <- as.numeric(as.character(theBed[,2]))
theBed[,3] <- as.numeric(as.character(theBed[,3]))
theBed[,5] <- as.numeric(as.character(theBed[,5]))

theSummits <- theBed[,2] + theBed[,5] - 1
theBed[,2] <- theSummits - as.numeric(args[2])
theBed[,3] <- theSummits + as.numeric(args[2])


theBed[theBed[,2]<1,2] <- 1

write.table(theBed, args[1], col.names=F, row.names=F, sep="\t", quote=F)
