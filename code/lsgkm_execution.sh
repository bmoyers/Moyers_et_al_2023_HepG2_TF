#!/bin/bash
#lsgkm_execution.sh

################################################################################
# The script below was used to produce gkmsvm models for use in downstream scripts,
#     specifically Running_gkmpredict.R .  The script below is tailored to
#     our local compute cluster's specific architecture, and will need to be
#     rewritten for your own compute environment's specific needs if you would
#     like to reproduce its results.  In particular, the various sourcings
#     will likely need to be removed, as well as the log_msg commands.  One will
#     need to provide or define their own TEMP_DIR.
#
# Expected to be in the user's PATH:
#     python 2.7.15
#     R 3.4.3
#     bedtools 2.28.0
#     lsgkm , see https://github.com/Dongwon-Lee/lsgkm
#     nullseq_generate.py , see https://beerlab.org/kmersvm/
#     Restrict_Peak_Size_summitCenter.R , provided in code.
#
# Required Environment Variables:
#     BED_FILE            - gzipped Bed file containing the appropriate IDR peaks for analysis.
#                             this will be restricted to the top 10k peaks (maximum) for
#                             speed and computational efficiency.
#     KMER_FILE           - A FASTA file containing all Kmers of interest.  Maximum Kmer and
#                             recommended Kmer size is 11mer. Provided as "All8Mers.txt"
#     HOT_SITES           - A bed-like file identifying all HOT sites, provided as "HOT_Sites.bed"
#     NULLSEQ_INDICES     - A directory containing nullseq indices for hg38.
#     GENOME_FILE         - A fasta for hg38.
#
# Optional Environment Variables:
#     OUTPUT_DIR          - Change default output directory
#     EXCLUSIONLIST           - An optional bedfile to remove
#
################################################################################

export SOURCE_DIR="/gpfs/gpfs2/sdi/etc"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

### Verify we are not running on the head node. ###
if [ -z "$LSB_JOBID" ]; then echo "Please run on a compute node. Exiting"; exit 1; fi

### Verify library directory exists. ###
if [ -z "$OUTPUT_DIR" ]; then OUTPUT_DIR=$(get_output_dir $LIBRARY);  fi
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p $OUTPUT_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $OUTPUT_DIR. Exiting"; exit 1; fi
fi


### Set up temp directory variable ###
if [ -z "$TEMP_DIR" ]; then TEMP_DIR=$(get_temp_dir); fi
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p $TEMP_DIR
    if [ $? -ne 0 ]; then echo "Could not create output dir: $TEMP_DIR. Exiting"; exit 1; fi
fi


log_msg info "Beginning lsgkm Analysis"
log_msg info "    BED:              $BED_FILE"
log_msg info "    OUTPUT_DIR:       $OUTPUT_DIR"

### Set up output files
BED_TEMP=$TEMP_DIR/Peaks_Temp.bed
BED_SORTED_TEMP1=$TEMP_DIR/Peaks_Temp_Sorted.bed
BED_SORTED_TEMP2=$TEMP_DIR/Peaks_Temp_Removing.bed
FASTA_TEMP=$TEMP_DIR/Peaks_Temp_top2k.fasta
NULL_MATCHED_BED=$TEMP_DIR/Null_matched.bed
NULL_MATCHED_FASTA=$TEMP_DIR/Null_matched.fasta
OUTPUT_BASE=$OUTPUT_DIR/gkmtrain_l8_k5_d2
OUTPUT_KMER=${OUTPUT_DIR}KmerEnrichments.txt


log_msg info "Unzipping bed file and restricting peaks..."

if file $BED_FILE | grep "gzip"
  then zcat $BED_FILE > $BED_TEMP
  else cp $BED_FILE $BED_TEMP
fi


CMD="sort -r -n -k 7 $BED_TEMP | cut -f 1,2,3,7,10 > $BED_SORTED_TEMP1"
run_cmd "$CMD" "$BED_SORTED_TEMP1"


grep -v "chrM" $BED_SORTED_TEMP1 | grep -v "random" | grep -v "EBV" | grep -v "chrUn" > $BED_SORTED_TEMP2
mv $BED_SORTED_TEMP2 $BED_SORTED_TEMP1

bedtools subtract -A -a $BED_SORTED_TEMP1 -b ${HOT_SITES} > $BED_SORTED_TEMP2
mv $BED_SORTED_TEMP2 $BED_SORTED_TEMP1


#If a exclusionList is provided, remove features which overlap with it.
if [[ -v EXCLUSIONLIST ]]
then
  bedtools subtract -A -a ${BED_SORTED_TEMP1} -b ${EXCLUSIONLIST} > $BED_SORTED_TEMP2
  mv $BED_SORTED_TEMP2 $BED_SORTED_TEMP1
else
  echo "ExclusionList not provided"
fi


#Restrict to the top 10k peaks
head -10000 $BED_SORTED_TEMP1 > $BED_SORTED_TEMP2
mv $BED_SORTED_TEMP2 $BED_SORTED_TEMP1


Rscript Restrict_Peak_Size_summitCenter.R $BED_SORTED_TEMP1 125


log_msg info "Generating Null Matched Sequences..."
CMD="python nullseq_generate.py -G -x 2 -m 1000 -r 1 -o $NULL_MATCHED_BED $BED_SORTED_TEMP1 hg38 ${NULLSEQ_INDICES}"
run_cmd "$CMD" "$NULL_MATCHED_BED"

log_msg info "Generating FASTA files from Bed files..."
CMD="fastaFromBed -fo $FASTA_TEMP -fi ${GENOME_FILE} -bed $BED_SORTED_TEMP1"
run_cmd "$CMD" "$FASTA_TEMP"

CMD="fastaFromBed -fo $NULL_MATCHED_FASTA -fi ${GENOME_FILE} -bed $NULL_MATCHED_BED"
run_cmd "$CMD" "$NULL_MATCHED_FASTA"



log_msg info "Running lsgkm..."
CMD="gkmtrain -T 16 -x 5 -r 1 -m 4000 -l 8 -k 5 -d 2 $FASTA_TEMP $NULL_MATCHED_FASTA $OUTPUT_BASE"
run_cmd "$CMD" "$OUTPUT_BASE.cvpred.txt"

CMD="gkmtrain -T 16 -m 4000 -l 8 -k 5 -d 2 $FASTA_TEMP $NULL_MATCHED_FASTA $OUTPUT_BASE"
run_cmd "$CMD" "$OUTPUT_BASE.model.txt"

log_msg info "Calculating Kmer enrichments..."
CMD="gkmpredict -T 16 $KMER_FILE $OUTPUT_BASE.model.txt $OUTPUT_KMER"

run_cmd "$CMD" "$OUTPUT_KMER"

log_msg info "Cleaning up files..."
CMD="rm $BED_TEMP $BED_TOP2K_TEMP1 $BED_TOP2K_TEMP2 $FASTA_TEMP $NULL_MATCHED_BED $NULL_MATCHED_FASTA"
run_cmd "$CMD"

### Finished
log_msg info "lsgkm_execution complete. Output is in $OUTPUT_DIR"
