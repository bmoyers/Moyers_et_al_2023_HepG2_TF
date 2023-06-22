#!/bin/bash
#Figure_4B.sh

################################################################################
#This script is used to produce Figure 4B.
#
#Requires that deepTools be in your PATH.
#This script takes as arguments:
# outDir: Path where files are to be written
# smaf_smaf: Path to smaf_smaf cobound regions, provided as sMAF_sMAF_cobound_peaks.bed
# smaf_cofa: Path to smaf_cofactor cobound regions, provided as sMAF_Cofactor_cobound_peaks.bed
# smaf_other: Path to smaf_smaf cobound regions, provided as sMAF_Other_cobound_peaks.bed
# other_other: Path to smaf_smaf cobound regions, provided as other_other_cobound_peaks.bed
# atacBW: Path to Bigwig of ATAC seq in HepG2, under ENCODE accession ENCFF262URW
################################################################################


###############################################################################
#Use deepTools to calculate ATAC-seq signal over the different peak sets and visualize
###############################################################################

outDir=${1}
smaf_smaf=${2}
smaf_cofa=${3}
smaf_other=${4}
other_other=${5}
atacBW=${6}

export logFile=${outDir}/deeptools_computeMatrix_logFile.txt

export outFile_smaf_atac=${outDir}/smaf_atac_matrix_all.mat.gz
computeMatrix reference-point  --referencePoint center -S ${atacBW} -R ${smaf_smaf} ${smaf_cofa} ${smaf_other} ${other_other} -a 500 -b 500 --outFileName ${outFile_smaf_atac}

export outFile_smaf_peak=${outDir}/Moyers_Figure4B.pdf
plotHeatmap -m ${outFile_smaf_atac} -out ${outFile_smaf_peak} --refPointLabel "Center" --regionsLabel s_s s_c s_o o_o --legendLocation center-right
