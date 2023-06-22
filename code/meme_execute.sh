#!/bin/bash
#meme_execute.sh

################################################################################
# The script below was used to produce meme outputs for use in downstream scripts,
#     The script below is tailored to our local compute cluster's specific
#     architecture, and will need to be rewritten for your own compute environment's
#     specific needs if you would like to reproduce its results.  In particular,
#     the various sourcings will likely need to be removed, as well as the log_msg
#     commands.  One will need to provide or define their own TEMP_DIR.
#
# Expected to be in the user's PATH:
#     meme 5.1.0
#     bedtools 2.28.0
#     python 2.7.15
#     R 3.4.3
#     Restrict_Peak_Size_summitCenter.R , provided in code.
#     Get_Flanking_Peak_Regions.R , provided in code
#     Motif_Control_Tests.R , provided in code.
#     MEME_Motif_Header.txt , provided in data, for construction of final files.
#
# Required Environment Variables:
#     BED_FILE            - gzipped Bed file containing the appropriate IDR peaks for analysis.
#                             this will be restricted to the top 10k peaks (maximum) for
#                             speed and computational efficiency.
#     BFILE               - A background file containing relevant background percentages for
#                             each nucleotide. provided as mammal_homo_sapiens_1000_199.na.bfile
#     GENOME_FILE         - A fasta for hg38.
#     HOT_SITES           - A bed-like file identifying all HOT sites, provided as "HOT_Sites.bed"
#     JASPAR_MOTIFS       - the JASPAR motif database in MEME motif format.
#     CISBP_MOTIFS        - the CISBP motif database in MEME motif format.
#
# Optional Environment Variables:
#     OUTPUT_DIR          - Change default output directory
#     BLACKLIST           - An optional bedfile to remove
#
#bsub -n 1 -R rusage[mem=6000] -We 20:00 -q c7normal -o /gpfs/gpfs1/home/bmoyers/ENCODE_lsgkm_meme_SMRTS/memeTest.bsubLog.txt -J meme "/gpfs/gpfs1/home/bmoyers/ENCODE_lsgkm_meme_SMRTS/Scripts/ENCODE_meme_pipeline_v2_2020March06.sh"
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


### Set up output files
BED_TEMP=$TEMP_DIR/Peaks_Temp.bed
BED_SORTED_TEMP1=$TEMP_DIR/Peaks_sorted_Temp1.bed
BED_TOP500_TEMP2=$TEMP_DIR/Peaks_top500_Temp2.bed
BED_SORTED_TEMP3=$TEMP_DIR/Peaks_Temp_removeProblems.bed


log_msg info "Unzipping bed file and restricting peaks..."

if file $BED_FILE | grep "gzip"
  then zcat $BED_FILE > $BED_TEMP
  else cp $BED_FILE $BED_TEMP
fi



CMD="sort -r -n -k 7 $BED_TEMP | cut -f 1,2,3,7 > $BED_SORTED_TEMP1"
run_cmd "$CMD" "$BED_SORTED_TEMP1"

grep -v "chrM" $BED_SORTED_TEMP1 | grep -v "random" | grep -v "EBV" | grep -v "chrUn" > $BED_SORTED_TEMP3
mv $BED_SORTED_TEMP3 $BED_SORTED_TEMP1


#If a Blacklist is provided, remove features which overlap with it.
if [[ -v BLACKLIST ]]
then
  bedtools subtract -A -a ${BED_SORTED_TEMP1} -b ${BLACKLIST} > $BED_SORTED_TEMP3
  mv $BED_SORTED_TEMP3 $BED_SORTED_TEMP1
else
  echo "Blacklist not provided"
fi

bedtools subtract -A -a $BED_SORTED_TEMP1 -b ${HOT_SITES} > $BED_SORTED_TEMP2
mv $BED_SORTED_TEMP2 $BED_SORTED_TEMP1


CMD="head -500 $BED_SORTED_TEMP1 > $BED_TOP500_TEMP2"
run_cmd "$CMD" "$BED_TOP500_TEMP2"

Rscript Restrict_Peak_Size.R $BED_TOP500_TEMP2 50




NUM_PEAKS=$(< "$BED_SORTED_TEMP1" wc -l)

if [ $NUM_PEAKS -ge 1000 ]; then

  FASTA_TOP500=${TEMP_DIR}/Peaks_for_motif_calling.fasta
  MOTIF_OUT=${OUTPUT_DIR}meme_chip_results

  BED_T1=${TEMP_DIR}/Peaks_T1.bed
  BED_C1=${TEMP_DIR}/Peaks_C1.bed
  BED_T2=${TEMP_DIR}/Peaks_T2.bed
  BED_C2=${TEMP_DIR}/Peaks_C2.bed

  FASTA_T1=${TEMP_DIR}/Peaks_T1.fasta
  FASTA_C1=${TEMP_DIR}/Peaks_C1.fasta
  FASTA_T2=${TEMP_DIR}/Peaks_T2.fasta
  FASTA_C2=${TEMP_DIR}/Peaks_C2.fasta

  FIMO_T1=${TEMP_DIR}/fimo_T1
  FIMO_C1=${TEMP_DIR}/fimo_C1
  FIMO_T2=${TEMP_DIR}/fimo_T2
  FIMO_C2=${TEMP_DIR}/fimo_C2

  log_msg info "Generating test sets, nullseqs, and FASTA files..."
  CMD="fastaFromBed -fo $FASTA_TOP500 -fi ${GENOME_FILE} -bed $BED_TOP500_TEMP2"
  run_cmd "$CMD" "$FASTA_TOP500"

  CMD="head -1000 $BED_SORTED_TEMP1 | tail -500 > $BED_T1"
  run_cmd "$CMD" "$BED_T1"
  Rscript Restrict_Peak_Size.R $BED_T1 150
  CMD="fastaFromBed -fo $FASTA_T1 -fi ${GENOME_FILE} -bed $BED_T1"
  run_cmd "$CMD" "$FASTA_T1"

  CMD="python /gpfs/gpfs1/home/bmoyers/Biotrain_2019/Scripts/nullseq_generate.py -x 1 -m 1000 -r 1 -o $BED_C1 $BED_TOP500_TEMP2 hg38 /gpfs/gpfs1/myerslab/reference/nullseq_indices/hg38_indices/"
  run_cmd "$CMD" "$BED_C1"
  Rscript Restrict_Peak_Size.R $BED_C1 150
  CMD="fastaFromBed -fo $FASTA_C1 -fi ${GENOME_FILE} -bed $BED_C1"
  run_cmd "$CMD" "$FASTA_C1"

  CMD="tail -n +501 $BED_SORTED_TEMP1 > $BED_T2"
  run_cmd "$CMD" "$BED_T2"
  Rscript Restrict_Peak_Size.R $BED_T2 150
  CMD="fastaFromBed -fo $FASTA_T2 -fi ${GENOME_FILE} -bed $BED_T2"
  run_cmd "$CMD" "$FASTA_T2"

  CMD="Rscript Get_Flanking_Peak_Regions.R $BED_T2 $BED_C2 300"
  run_cmd "$CMD" "$BED_C2"
  CMD="fastaFromBed -fo $FASTA_C2 -fi ${GENOME_FILE} -bed $BED_C2"
  run_cmd "$CMD" "$FASTA_C2"



  log_msg info "Running MEME-ChIP..."
  CMD="meme-chip -oc $MOTIF_OUT -noecho -bfile $BFILE -dna -meme-mod zoops -meme-nmotifs 5 -meme-minw 6 -meme-maxw 50 -spamo-skip -fimo-skip $FASTA_TOP500"
  MEME_OUTPUT_FILE="${MOTIF_OUT}/meme_out/meme.txt"
  run_cmd "$CMD" "$MEME_OUTPUT_FILE"

  log_msg info "Running FIMO controls..."

  fimo --oc $FIMO_T1 $MEME_OUTPUT_FILE $FASTA_T1
  fimo --oc $FIMO_C1 $MEME_OUTPUT_FILE $FASTA_C1
  fimo --oc $FIMO_T2 $MEME_OUTPUT_FILE $FASTA_T2
  fimo --oc $FIMO_C2 $MEME_OUTPUT_FILE $FASTA_C2

  PASSED_MOTIFS=${OUTPUT_DIR}Passed_Motifs
  mkdir $PASSED_MOTIFS
  PARTIAL_MOTIFS=${OUTPUT_DIR}Partial_Motifs
  mkdir $PARTIAL_MOTIFS
  FAILED_MOTIFS=${OUTPUT_DIR}Failed_Motifs
  mkdir $FAILED_MOTIFS

  Rscript Motif_Control_Tests.R $MEME_OUTPUT_FILE $FIMO_T1/fimo.tsv $FIMO_C1/fimo.tsv $FIMO_T2/fimo.tsv $FIMO_C1/fimo.tsv $PASSED_MOTIFS $PARTIAL_MOTIFS $FAILED_MOTIFS $BED_FILE $BED_C1 $BED_T2 $BED_C2 $MOTIF_OUT/meme_out ${MOTIF_OUT}/centrimo_out/centrimo.tsv

  ALL_PASSED_MOTIFS=${OUTPUT_DIR}All_Passed_Motifs.txt
  PARTIAL_PASSED_MOTIFS=${OUTPUT_DIR}Partial_Passed_Motifs.txt
  ALL_MOTIFS=${OUTPUT_DIR}All_Motifs.txt
  cat MEME_Motif_Header.txt $PASSED_MOTIFS/*.txt > $ALL_PASSED_MOTIFS
  cat MEME_Motif_Header.txt $PARTIAL_MOTIFS/*.txt > $PARTIAL_PASSED_MOTIFS
  cat MEME_Motif_Header.txt $PASSED_MOTIFS/*.txt $PARTIAL_MOTIFS/*.txt  $FAILED_MOTIFS/*.txt > $ALL_MOTIFS


  if  cat ${ALL_MOTIFS} | grep "MOTIF" ; then

    log_msg info "Running TOMTOM and FIMO..."
    CMD="tomtom -oc ${OUTPUT_DIR}JASPAR_TOMTOM/  $ALL_MOTIFS ${JASPAR_MOTIFS}"
    run_cmd "$CMD" "${OUTPUT_DIR}JASPAR_TOMTOM/tomtom.tsv"

    CMD="tomtom -oc ${OUTPUT_DIR}CISBP_TOMTOM/ $ALL_MOTIFS ${CISBP_MOTIFS}"
    run_cmd "$CMD" "${OUTPUT_DIR}CISBP_TOMTOM/tomtom.tsv"

    FASTA_FIMO=${TEMP_DIR}/AllPeaks.fasta
    CMD="fastaFromBed -fo $FASTA_FIMO -fi ${GENOME_FILE} -bed $BED_TEMP"
    run_cmd "$CMD" "$FASTA_FIMO"

    log_msg info "Running Fimo..."
    CMD="fimo --oc ${OUTPUT_DIR}FIMO_Peaks/ ${ALL_MOTIFS} ${FASTA_FIMO}"
    run_cmd "$CMD" "${OUTPUT_DIR}FIMO_Peaks/fimo.tsv"


    log_msg info "ENCODE_meme_pipeline complete. Output is in $MOTIF_OUT"


  else log_msg info "No interesting motifs.  TomTom and FIMO not run."

  fi

fi



if [ $NUM_PEAKS -lt 1000 ]; then
  log_msg info "Fewer than 1000 peaks!  Meme analyses were not run."
fi

log_msg info "Cleaning up files..."
CMD="rm -rf ${TEMP_DIR}"
run_cmd "$CMD"
