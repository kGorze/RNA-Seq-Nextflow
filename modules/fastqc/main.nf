#!/usr/bin/env nextflow

/*
 * Enhanced FastQC module for RNA-Seq pipeline
 * Includes adapter detection and comprehensive QC metrics
 */

// Import adapter file creation process
include { CREATE_ADAPTER_FILE } from './adapters'

workflow FASTQC_WORKFLOW {
    take:
    reads_ch  // Channel with sample_id and reads

    main:
    // Create adapter file
    CREATE_ADAPTER_FILE()
    adapter_file = CREATE_ADAPTER_FILE.out.adapter_file
    
    // Run FastQC
    FASTQC(reads_ch, adapter_file)
    
    // Output channels
    fastqc_zip = FASTQC.out.zip
    fastqc_html = FASTQC.out.html
    
    emit:
    zip = fastqc_zip
    html = fastqc_html
}

process FASTQC {
    tag "$sample_id"
    label 'process_low'
    
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path adapter_file
    
    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html
    
    script:
    """
    # Create directory for adapter file
    mkdir -p \$HOME/.fastqc/
    cp ${adapter_file} \$HOME/.fastqc/adapters.fa
    
    # Run FastQC with adapter detection
    fastqc -q -t ${task.cpus} \\
        --adapters \$HOME/.fastqc/adapters.fa \\
        --extract \\
        --nogroup \\
        ${reads}
    """
}
