/*
 * pipeline input parameters
 */
params.reads = "$projectDir/raw_reads/*_{1,2}.fastq.gz"
params.transcriptome_file = "$projectDir/ref/dmel-all-chromosome-r6.46.fasta"
params.annotationfile_file = "$projectDir/ref/dmel-all-r6.46.gtf"
params.outdir = "results"
log.info """\
    T R A N S E Q F L O W - P E
    ===========================
    Author       : Saifeldeen M. Ibrahim
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

process FASTQC {
    tag "FASTQC: $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

process TRIM {
    tag "Trimmomatic: $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed_reads/Paired/${sample_id}_1_paired.fastq"), path("trimmed_reads/Paired/${sample_id}_2_paired.fastq")

    script:
    """
    mkdir -p trimmed_reads/Paired trimmed_reads/Unpaired
    TrimmomaticPE -phred33 ${reads[0]} ${reads[1]}\
    trimmed_reads/Paired/${sample_id}_1_paired.fastq trimmed_reads/Unpaired/${sample_id}_1_unpaired.fastq\
    trimmed_reads/Paired/${sample_id}_2_paired.fastq trimmed_reads/Unpaired/${sample_id}_2_unpaired.fastq\
    SLIDINGWINDOW:4:25 MINLEN:36
    """
}

process INDEX {
    tag "Reference Indexing"
    publishDir params.outdir, mode:'copy'

    input:
    path 'refrence'
    path 'GTF'

    output:
    path "ref_index"

    script:
    """
    STAR --runThreadN 2\
     --runMode genomeGenerate\
      --genomeDir ${"ref_index"}\
       --genomeFastaFiles ${'refrence'}\
        --sjdbGTFfile ${'GTF'}\
         --genomeSAindexNbases 14
    """
}

process MAPPING {
    tag "STAR-Alignment: $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path 'index'
    path 'GTF'
    
    output:
    tuple val(sample_id), path("STAR-Aligned/${sample_id}.sam")

    script:
    """
    mkdir -p STAR-Aligned
    STAR --runThreadN 2\
     --genomeDir ${'index'}\
      --sjdbGTFfile ${'GTF'}\
       --readFilesIn ${read1} ${read2}\
        #if your reads are in gz format
         #--readFilesCommand zcat\
          --outSAMtype BAM SortedByCoordinate\
           --limitBAMsortRAM 3000000000\
            --quantMode GeneCounts
    mv Aligned.out.sam ${sample_id}.sam
    mv *.sam STAR-Aligned
    mv *.out STAR-Aligned
    mv *.tab STAR-Aligned

    """
}

process BAM {
    tag "Convert SAM to BAM: $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(sam)

    
    output:
    tuple val(sample_id), path("BAM-files/${sample_id}.bam")

    script:
    """
    mkdir -p BAM-files
    # convert sam to bam by view and sort by samtools
    samtools view -bS ${sam} > ${sample_id}.bam
    # sort the bam file
    samtools sort ${sample_id}.bam -o ${sample_id}_sorted.bam
    # index the bam file
    samtools index ${sample_id}_sorted.bam
    # remove sam and bam files
    mv ${sample_id}_sorted.bam BAM-files/${sample_id}.bam
    """
}

process FeatureCount {
    tag "FeatureCount: $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path 'GTF'
    tuple val(sample_id), path(bam)

    
    output:
    tuple val(sample_id), path("Counts/${sample_id}.txt")

    script:
    """
    mkdir -p Counts
    featureCounts -p -a ${GTF} -o Counts/${sample_id}_out.txt ${bam}
    awk '{print \$1, \$NF}' Counts/${sample_id}_out.txt > Counts/${sample_id}.txt
    """
}


workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads)
    fastqc_ch = FASTQC(read_pairs_ch)
    multiqc_ch = MULTIQC(fastqc_ch.collect())
    trimmomatic_ch = TRIM(read_pairs_ch)
    index_ch = INDEX(params.transcriptome_file, params.annotationfile_file)
    mapping_ch = MAPPING(trimmomatic_ch, index_ch, params.annotationfile_file)
    bamfile_ch = BAM(mapping_ch)
    counts_ch = FeatureCount(params.annotationfile_file, bamfile_ch)
}
workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Explore the results at --> $params.outdir\n" : "Oops .. something went wrong" )
}
