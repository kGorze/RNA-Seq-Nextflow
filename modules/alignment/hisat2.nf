#!/usr/bin/env nextflow

/*
 * Enhanced HISAT2 alignment workflow for RNA-Seq pipeline
 * Includes genome indexing, alignment, and BAM processing
 */

workflow HISAT2 {
    take:
    reads_ch     // Channel with sample_id and reads
    genome_file  // Reference genome FASTA
    gtf_file     // Gene annotation GTF
    
    main:
    // Extract splice sites and exons from GTF
    EXTRACT_SPLICE_SITES(gtf_file)
    EXTRACT_EXONS(gtf_file)
    
    // Generate HISAT2 index
    HISAT2_INDEX(
        genome_file, 
        EXTRACT_SPLICE_SITES.out.splice_sites, 
        EXTRACT_EXONS.out.exons
    )
    
    // Align reads with HISAT2
    HISAT2_ALIGN(
        reads_ch, 
        HISAT2_INDEX.out.index, 
        EXTRACT_SPLICE_SITES.out.splice_sites
    )
    
    // Sort and index BAM files
    SAMTOOLS_SORT_INDEX(HISAT2_ALIGN.out.bam)
    
    // Generate alignment statistics
    SAMTOOLS_STATS(SAMTOOLS_SORT_INDEX.out.sorted_bam)
    
    emit:
    bam = SAMTOOLS_SORT_INDEX.out.sorted_bam
    bai = SAMTOOLS_SORT_INDEX.out.bam_index
    stats = SAMTOOLS_STATS.out.stats
    logs = HISAT2_ALIGN.out.log
}

process EXTRACT_SPLICE_SITES {
    label 'process_low'
    
    container 'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3'
    
    input:
    path gtf
    
    output:
    path "splicesites.txt", emit: splice_sites
    
    script:
    """
    hisat2_extract_splice_sites.py ${gtf} > splicesites.txt
    """
}

process EXTRACT_EXONS {
    label 'process_low'
    
    container 'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3'
    
    input:
    path gtf
    
    output:
    path "exons.txt", emit: exons
    
    script:
    """
    hisat2_extract_exons.py ${gtf} > exons.txt
    """
}

process HISAT2_INDEX {
    label 'process_high'
    
    container 'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3'
    
    publishDir "${params.outdir}/hisat2/index", mode: 'copy'
    
    input:
    path genome
    path splice_sites
    path exons
    
    output:
    path "hisat2_index", emit: index
    
    script:
    """
    mkdir -p hisat2_index
    
    hisat2-build \\
        -p ${task.cpus} \\
        --ss ${splice_sites} \\
        --exon ${exons} \\
        ${genome} \\
        hisat2_index/genome
    
    # Move all index files to the index directory
    mv genome.*.ht2 hisat2_index/ || true
    """
}

process HISAT2_ALIGN {
    tag "$sample_id"
    label 'process_high'
    
    container 'quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3'
    
    publishDir "${params.outdir}/hisat2/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path index
    path splice_sites
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}.hisat2.log", emit: log
    
    script:
    def prefix = "${sample_id}"
    def strandedness = ""
    if (params.strandedness == 'forward') {
        strandedness = "--rna-strandness FR"
    } else if (params.strandedness == 'reverse') {
        strandedness = "--rna-strandness RF"
    }
    
    def read_args = ""
    if (params.paired_end) {
        read_args = "-1 ${reads[0]} -2 ${reads[1]}"
    } else {
        read_args = "-U ${reads}"
    }
    
    """
    hisat2 \\
        -p ${task.cpus} \\
        -x ${index}/genome \\
        ${read_args} \\
        --known-splicesite-infile ${splice_sites} \\
        --novel-splicesite-outfile ${prefix}.novel_splicesites.txt \\
        --summary-file ${prefix}.hisat2.log \\
        --dta \\
        ${strandedness} \\
        --met-file ${prefix}.hisat2.metrics.txt \\
        | samtools view -bS - > ${prefix}.bam
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/samtools:1.13--h8c37831_0'
    
    publishDir "${params.outdir}/hisat2/sorted", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam
    path "${sample_id}.sorted.bam.bai", emit: bam_index
    
    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}

process SAMTOOLS_STATS {
    tag "$sample_id"
    label 'process_low'
    
    container 'quay.io/biocontainers/samtools:1.13--h8c37831_0'
    
    publishDir "${params.outdir}/hisat2/stats", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path "${sample_id}.stats.txt", emit: stats
    path "${sample_id}.flagstat.txt", emit: flagstat
    path "${sample_id}.idxstats.txt", emit: idxstats
    
    script:
    """
    samtools stats ${bam} > ${sample_id}.stats.txt
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt
    """
}
