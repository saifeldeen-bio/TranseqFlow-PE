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

#### PS: you have to make sure you have installed NextFlow and all of these tools on your environment before using this pipline.

# Input Parameters

The pipeline accepts the following input parameters:
  
  `reads`: Path to paired-end reads in .fastq.gz format. The pipeline expects the reads to be organized as forward and reverse reads for each sample.
  
  `transcriptome_file`: Path to the reference transcriptome file in FASTA format.
  
  `annotationfile_file`: Path to the reference genome annotation file in GTF format.
  
  `outdir`: Directory where the pipeline will store the output files.

#### Please make sure that your fastqs are in .fastq.gz format and all of reads are placed in `raw_reads` directory, and your reference and gtf files are stored in `ref` directory. 

# Usage

```bash
nextflow run transseqflowPE --reads '/path/to/raw_reads/*_{1,2}.fastq.gz' --transcriptome_file '/path/to/ref/reference.fasta' --annotationfile_file '/path/to/ref/annotation.gtf'
```

# Outputs

The results are saved in the specified output directory (outdir) and include:

  - FASTQC quality reports

  - A MultiQC summary report

  - Trimmed reads

  - BAM files (sorted and indexed)

  - Gene count tables for each sample
