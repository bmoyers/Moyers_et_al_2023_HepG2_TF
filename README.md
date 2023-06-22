# Moyers et al 2023 HepG2 TF

Welcome to the code and data repository for the forthcoming paper, "Characterization of human transcription actor function and patterns of gene regulation in HepG2 cells".

This repository has a simple structure, with the "code" folder containing relevant code for understanding and replicating key analyses in the paper.  Data used in the execution for these scripts is noted at the header of each one.  Much of this data is provided in the "data" folder.  However, some data is not provided as it is available through other means, such as ENCODE accessions, which can be accessed on https://www.encodeproject.org/, or motif repositories available at https://jaspar.genereg.net/ or http://cisbp.ccbr.utoronto.ca/.

Find below a description for each item provided in the data folder.  If you have any questions,
please feel free to contact me.

- Agarwal2023: Directory containing processed MPRA data from [Massively parallel characterization of transcriptional regulatory elements in three diverse human cell types](https://www.biorxiv.org/content/10.1101/2023.03.05.531189v1.full)
- All8Mers.txt: a FASTA file containing all possible 8mers, used for construction of gkmsvm models.
- All_Merged_Peaks.bed: A bed-like file containing all peaks, sorted and merged, from all "preferred" experiments. Column 4 contains the total number of DAPs covering a given entry, and Column 5 contains the identity of each DAP contained within a region.
- Binding_Expression_Matrix.txt: Pre-computed matrix which relates, for each gene, whether or not a given TF has a binding within +/- 500 bp of a gene's TSS.
- Binding_Expression_Matrix_anyCGI.txt: As "Binding_Expression_Matrix", but only for genes with a CpG Island in the promoter.
- Binding_Expression_Matrix_nonCGI.txt: As "Binding_Expression_Matrix", but only for genes without a CpG Island in the promoter.
- Element_representation_for_Rep1.txt, Rep2, Rep3: Processed counts data and signal in each replicate for the number of reads in RNA and DNA for each element in our lentiviral MPRA. Raw and processed data files can be found at [GSE235360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235360)
- Experiment_Set: Directory containing, for each ChIP-seq dataset in HepG2, the processed narrowPeak files. For experiments released on the [ENCODE](https://www.encodeproject.org/) browser, accession numbers are provided.  For those listed as "Local", raw and processed files can be found at [GSE235477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235477)
- Experiment_Set_K562: As Experiment_Set, but for ChIP-seq datasets in K562 which were released on ENCODE. These datasets were only used for cross-cell-type validation of Binding-expression models.
- HOT_Sites.bed: As All_Merged_Peaks.bed, but containing only regions classified as HOT sites.
- MEME_Motif_Header.txt: A simple meme header used for production of final files by the meme_execute.sh script.
- binding_expr_models_results.rds: Data output produced by the binding_expr_modeling.R script.
- finding_best_model_size.rds: Data output produced by the find_model_size.R script.
- gkmsvm_models.tar.gz: A gzipped directory of gkmsvm models produced by the lsgkm_execution.sh script.
- hg38.chrom.sizes: A tab-delimited file denoting the chromosome sizes for hg38.
- mammal_homo_sapiens_1000_199.na.bfile: A file used by the meme_execute.sh script which contains nucleotide frequencies.
- meme_passed_motifs: A directory containing the passed motifs for each ChIP-seq dataset produced by the meme_execute.sh script.
- non_cCRE_allFactors_regions_with_various_annotations.txt: a tab-delimited file containing non-cCRE regions which were bound by 2 or more DAPs in our datasets, annotated with various genomic features.
- other_other_cobound_peaks.bed: A bed-like file containing a subsampling of regions which were bound by 2 or more factors (excluding sMAF factors and their cofactors) with peak centers within 50bp of one another.
- proms_v103_annotated_GR_like.tsv: A tab delimited file which annotates each gene as having or not having CpG Islands in their promoters.  
- refseq_genes_unique_TSS.bed: refseq TSSes.
- refseq_genes_unique_TSS_1000.bed: refseq TSSes +/- 500bp.
- sMAF_Cofactor_cobound_peaks.bed: A bed-like file containing regions which were bound by a sMAF and a noted sMAF cofactor with peak centers within 50bp of one another.
- sMAF_Other_cobound_peaks.bed: A bed-like file containing regions which were bound by a sMAF and another factor (excluding sMAFs and noted sMAF cofactors) with peak centers within 50bp of one another.
- sMAF_sMAF_cobound_peaks.bed: A bed-like file containing regions which were bound two sMAF factors with peak centers within 50bp of one another.
- unbound_regions_cCREs_in_ATAC.bed: A bed file containing unbound cCRE regions which were found to be within ATAC-seq peaks.  
