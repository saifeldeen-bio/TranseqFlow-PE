# TranseqFlow-PE
Transcriptome upstream analysis pipeline with Nextflow for paired-end reads. It performs quality check with FASTQC and MultiQC, trimming with Trimmomatic, alignment with STAR, and quantification of RNA-seq data with FeatureCount, producing clean and reproducible results. 

## The Processs and the inegrated tools
  `FASTQC - fastqc`: Performs quality control checks on raw sequencing data to generate quality reports.

  `MULTIQC - multiqc`: Aggregates results from FASTQC into a single report for easier analysis.

  `TRIM - trimmpmaticPE`: Uses Trimmomatic to remove low-quality bases and adapter sequences from paired-end reads, generating trimmed and cleaned reads.

  `INDEX - STAR`: Creates a genome index using the reference transcriptome and annotation file with STAR to prepare for the alignment step.

  `MAPPING - STAR`: Aligns the trimmed reads to the reference genome using STAR, producing SAM files for further processing.

  `BAM - samtools`: Converts SAM files to BAM format, sorts and indexes them for downstream analysis using samtools.

  `FeatureCount - featureCounts`: Quantifies gene expression by counting the number of reads aligned to each gene using featureCounts.

#### PS: you have to make sure you have installed all of these tools on your environment before using this pipline.

