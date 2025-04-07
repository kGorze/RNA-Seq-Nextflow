#!/usr/bin/env nextflow

/*
 * Enhanced STAR alignment workflow for RNA-Seq pipeline
 * Includes genome indexing, alignment, and BAM processing
 */

workflow STAR_WORKFLOW {
    take:
    reads_ch     // Channel with sample_id and reads
    genome_file  // Reference genome FASTA
    gtf_file     // Gene annotation GTF
    
    main:
    // Generate STAR index
    STAR_INDEX(genome_file, gtf_file)
    
    // Align reads with STAR
    STAR_ALIGN(reads_ch, STAR_INDEX.out.index, gtf_file)
    
    // Sort and index BAM files
    SAMTOOLS_SORT_INDEX(STAR_ALIGN.out.bam)
    
    // Generate alignment statistics
    SAMTOOLS_STATS(SAMTOOLS_SORT_INDEX.out.sorted_bam)
    
    emit:
    bam = SAMTOOLS_SORT_INDEX.out.sorted_bam
    bai = SAMTOOLS_SORT_INDEX.out.bam_index
    stats = SAMTOOLS_STATS.out.stats
    logs = STAR_ALIGN.out.log
}

process STAR_INDEX {
    label 'process_high'
    
    container 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'
    
    publishDir "${params.outdir}/star/index", mode: 'copy'
    
    input:
    path genome
    path gtf
    
    output:
    path "star_index", emit: index
    
    script:
    """
    mkdir -p star_index
    
    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100 \\
        --genomeSAindexNbases 11
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    label 'process_high'
    
    container 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'
    
    publishDir "${params.outdir}/star/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path index
    path gtf
    
    output:
    tuple val(sample_id), path("${sample_id}.Aligned.out.bam"), emit: bam
    path "*.Log.*", emit: log
    path "*.SJ.out.tab", emit: splice_junctions
    path "*.ReadsPerGene.out.tab", optional: true, emit: gene_counts
    
    script:
    def prefix = "${sample_id}"
    def read_args = params.paired_end ? "--readFilesIn ${reads[0]} ${reads[1]}" : "--readFilesIn ${reads}"
    def strandedness = ""
    if (params.strandedness == 'forward') {
        strandedness = "--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical"
    } else if (params.strandedness == 'reverse') {
        strandedness = "--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical"
    }
    
    """
    STAR --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesCommand zcat \\
        ${read_args} \\
        --outFileNamePrefix ${prefix}. \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes Standard \\
        --outSAMunmapped Within \\
        --outFilterType BySJout \\
        --outFilterMultimapNmax 20 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --sjdbGTFfile ${gtf} \\
        --quantMode GeneCounts \\
        ${strandedness}
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/samtools:1.13--h8c37831_0'
    
    publishDir "${params.outdir}/star/sorted", mode: 'copy'
    
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
    
    publishDir "${params.outdir}/star/stats", mode: 'copy'
    
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
