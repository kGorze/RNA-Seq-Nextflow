#!/usr/bin/env nextflow

/*
 * Trimmomatic module for RNA-Seq pipeline
 */

process TRIMMOMATIC {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
    
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: log
    
    script:
    def prefix = "${sample_id}"
    if (params.paired_end) {
        """
        trimmomatic PE -threads ${task.cpus} \\
            ${reads[0]} ${reads[1]} \\
            ${prefix}_1_trimmed.fastq.gz ${prefix}_1_unpaired.fastq.gz \\
            ${prefix}_2_trimmed.fastq.gz ${prefix}_2_unpaired.fastq.gz \\
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \\
            LEADING:3 TRAILING:3 MINLEN:36 \\
            2> ${prefix}_trimming_report.txt
        """
    } else {
        """
        trimmomatic SE -threads ${task.cpus} \\
            ${reads} \\
            ${prefix}_trimmed.fastq.gz \\
            ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \\
            LEADING:3 TRAILING:3 MINLEN:36 \\
            2> ${prefix}_trimming_report.txt
        """
    }
}
