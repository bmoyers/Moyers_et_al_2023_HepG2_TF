#!/usr/bin/env R
#Motif_Control_Tests.R


################################################################################
#This script takes as arguments several fimo and centrimo results produced during
#the running of meme_execute.sh , as well as directories to write various results to.
#It performs tests regarding motif count in various peak sets as well as
#centrality metrics to determine if a motif is passed, partially passed, or failed.
################################################################################

################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################

extract_motifs <- function(motifs_File, bed_file) {
  thisName <- strsplit(bed_file, split="/")[[1]]
  thisName <- thisName[length(thisName)]
  thisName <- gsub(".bed", "", thisName)
  thisName <- gsub(".gz", "", thisName)

  theseLines <- readLines(motifs_File)
  theseLines <- theseLines[-grep("^chr", theseLines)]

  theStarts <- grep("^letter-probability matrix", theseLines)
  theEnds <- grep("regular expression$", theseLines)

  allMotifs <- list()
  for (i in 1:length(theStarts)) {
    motifName <- paste("MOTIF ", thisName, "_motif_", i, " ", "MOTIF-", i, sep="")
    thisMotif <- theseLines[theStarts[i]:(theEnds[i]-4)]

    thisMotif <- c(motifName, thisMotif, "", "")
    allMotifs[[length(allMotifs)+1]] <- thisMotif
  }
  names(allMotifs) <- paste("MEME-", 1:length(theStarts), sep="")
  return(allMotifs)
}

extract_evals <- function(motifs_File) {
  theseLines <- readLines(motifs_File)
  theseLines <- theseLines[-grep("^chr", theseLines)]
  theseLines <- theseLines[grep("MOTIF", theseLines)]
  theseLines <- theseLines[grep("E-value", theseLines)]

  theEvals <- c()
  for (i in 1:length(theseLines)) {
    thisEval <- strsplit(theseLines[i], split=" ")[[1]]
    thisEval <- thisEval[length(thisEval)]
    thisEval <- as.numeric(as.character(thisEval))
    theEvals <- c(theEvals, thisEval)
  }
  return(theEvals)
}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)

motifs_File <- args[1]
fimo_t1 <- args[2]
fimo_c1 <- args[3]
fimo_t2 <- args[4]
fimo_c2 <- args[5]
success_dir <- args[6]
partial_dir <- args[7]
failed_dir <- args[8]
bed_file <- args[9]
bed_c1 <- args[10]
bed_t2 <- args[11]
bed_c2 <- args[12]
meme_out <- args[13]
centrimo <- args[14]



allMotifs <- extract_motifs(motifs_File, bed_file)
theEvals <- extract_evals(motifs_File)

fimo_t1 <- read.table(fimo_t1, header=T, sep="\t", stringsAsFactors=F)
fimo_c1 <- read.table(fimo_c1, header=T, sep="\t", stringsAsFactors=F)
fimo_t2 <- read.table(fimo_t2, header=T, sep="\t", stringsAsFactors=F)
fimo_c2 <- read.table(fimo_c2, header=T, sep="\t", stringsAsFactors=F)

bed_c1 <- read.table(bed_c1, header=F, sep="\t", stringsAsFactors=F)
bed_c1_names <- paste(bed_c1[,1], ":", bed_c1[,2], "-", bed_c1[,3], sep="")
bed_t2 <- read.table(bed_t2, header=F, sep="\t", stringsAsFactors=F)
bed_c2 <- read.table(bed_c2, header=F, sep="\t", stringsAsFactors=F)

logo_files <- list.files(meme_out, pattern="logo", full.names=T)

centrimo <- read.table(centrimo, header=T, sep="\t", stringsAsFactors=F)



motif_names <- names(allMotifs)
for (i in 1:length(motif_names)) {
  this_fimo_t1 <- fimo_t1[fimo_t1[,"motif_alt_id"]==motif_names[i],]
  T1 <- length(unique(as.character(this_fimo_t1[,"sequence_name"])))

  this_fimo_c1 <- fimo_c1[fimo_c1[,"motif_alt_id"]==motif_names[i],]
  these_C1_Nums <- c()
  for (j in 1:100) {
    theseSeqs <- sample(bed_c1_names, size=500, replace=T)
    these_C1_Nums <- c(these_C1_Nums, length(theseSeqs[theseSeqs%in%as.character(this_fimo_c1[,"sequence_name"])]))
  }
  mean_C1_nums <- mean(these_C1_Nums)
  sd_C1_nums <- sd(these_C1_Nums)
  zscore_T1 <- (T1-mean_C1_nums)/sd_C1_nums
  pval_T1 <- pnorm(-zscore_T1)
  qval_T1 <- p.adjust(pval_T1, method="fdr", n=5)


  this_fimo_t2 <- fimo_t2[fimo_t2[,"motif_alt_id"]==motif_names[i],]
  T2 <- length(unique(as.character(this_fimo_t2[,"sequence_name"]))) / nrow(bed_t2)
  this_fimo_c2 <- fimo_c2[fimo_c2[,"motif_alt_id"]==motif_names[i],]
  C2 <- length(unique(as.character(this_fimo_c2[,"sequence_name"]))) / nrow(bed_c2)


  T1_success <- FALSE
  if(qval_T1<0.00001) { T1_success <- TRUE }

  T2_success <- FALSE
  if((T2 > 0.10) && (T2/C2 >= 1.25)) {T2_success <- TRUE}

  T3_success <- FALSE
  if(motif_names[i]%in%centrimo[,"motif_alt_id"]) {
    this_eval <- as.numeric(as.character(centrimo[which(centrimo[,"motif_alt_id"]==motif_names[i]),"E.value"]))
    if(this_eval<0.05) { T3_success <- TRUE }
  }

  T4_success <- FALSE
  if(theEvals[i]<0.05) { T4_success <- TRUE }


  this_Motif_Name <- strsplit(allMotifs[[i]][1], split=" ")[[1]]
  this_Motif_Name <- this_Motif_Name[this_Motif_Name!=""]
  this_Motif_Name <- this_Motif_Name[2]


  FactorBook_Success <- FALSE
  if(T1_success && T2_success && T4_success) { FactorBook_Success <- TRUE }

  HudsonAlpha_Success <- FALSE
  if(T1_success && T3_success && T4_success) { HudsonAlpha_Success <- TRUE }


  if(FactorBook_Success && HudsonAlpha_Success ) {
    thisOutput <- paste(success_dir, "/", this_Motif_Name, ".txt", sep="")
    writeLines(allMotifs[[i]], thisOutput, sep="\n")
    these_logos <- logo_files[grep(paste(as.character(i), ".png", sep=""), logo_files)]
    for (j in 1:length(these_logos)) {
      logo_name <- strsplit(these_logos[j], split="/")[[1]]
      logo_name <- logo_name[length(logo_name)]
      theCommand <- paste("cp ", these_logos[j], " ", success_dir, "/", this_Motif_Name, "_", logo_name, sep="")
      system(theCommand)
    }
  }

  if(FactorBook_Success && !HudsonAlpha_Success) {
    thisOutput <- paste(partial_dir, "/", this_Motif_Name, "_Failed_HA_v7_test.txt", sep="")
    writeLines(allMotifs[[i]], thisOutput, sep="\n")
    these_logos <- logo_files[grep(paste(as.character(i), ".png", sep=""), logo_files)]
    for (j in 1:length(these_logos)) {
      logo_name <- strsplit(these_logos[j], split="/")[[1]]
      logo_name <- logo_name[length(logo_name)]
      theCommand <- paste("cp ", these_logos[j], " ", partial_dir, "/", this_Motif_Name, "_", logo_name, sep="")
      system(theCommand)
    }
  }

  if(!FactorBook_Success && HudsonAlpha_Success) {
    thisOutput <- paste(partial_dir, "/", this_Motif_Name, "_Failed_FactBook_v2_test.txt", sep="")
    writeLines(allMotifs[[i]], thisOutput, sep="\n")
    these_logos <- logo_files[grep(paste(as.character(i), ".png", sep=""), logo_files)]
    for (j in 1:length(these_logos)) {
      logo_name <- strsplit(these_logos[j], split="/")[[1]]
      logo_name <- logo_name[length(logo_name)]
      theCommand <- paste("cp ", these_logos[j], " ", partial_dir, "/", this_Motif_Name, "_", logo_name, sep="")
      system(theCommand)
    }
  }

  if(!FactorBook_Success && !HudsonAlpha_Success ) {
    thisOutput <- paste(failed_dir, "/", this_Motif_Name, ".txt", sep="")
    writeLines(allMotifs[[i]], thisOutput, sep="\n")
    these_logos <- logo_files[grep(paste(as.character(i), ".png", sep=""), logo_files)]
    for (j in 1:length(these_logos)) {
      logo_name <- strsplit(these_logos[j], split="/")[[1]]
      logo_name <- logo_name[length(logo_name)]
      theCommand <- paste("cp ", these_logos[j], " ", failed_dir, "/", this_Motif_Name, "_", logo_name, sep="")
      system(theCommand)
    }
  }

}
