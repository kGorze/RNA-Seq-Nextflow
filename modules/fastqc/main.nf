#!/usr/bin/env nextflow

/*
 * Enhanced FastQC module for RNA-Seq pipeline
 * Includes adapter detection and comprehensive QC metrics
 */

workflow FASTQC {
    take:
    reads_ch  // Channel with sample_id and reads

    main:
    // Create adapter file
    CREATE_ADAPTER_FILE()
    adapter_file = CREATE_ADAPTER_FILE.out.adapter_file
    
    // Run FastQC
    RUN_FASTQC(reads_ch, adapter_file)
    
    // Output channels
    fastqc_zip = RUN_FASTQC.out.zip
    fastqc_html = RUN_FASTQC.out.html
    
    emit:
    zip = fastqc_zip
    html = fastqc_html
}

process CREATE_ADAPTER_FILE {
    label 'process_low'
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    output:
    path "illumina_adapters.txt", emit: adapter_file
    
    script:
    """
    cat > illumina_adapters.txt << EOL
TruSeq_Universal_Adapter\tAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
TruSeq_Adapter_Index_1\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
TruSeq_Adapter_Index_2\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
TruSeq_Adapter_Index_3\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
TruSeq_Adapter_Index_4\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
TruSeq_Adapter_Index_5\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGAATCTCGTATGCCGTCTTCTGCTTG
TruSeq_Adapter_Index_6\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
TruSeq_Adapter_Index_7\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
TruSeq_Adapter_Index_8\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
EOL
    """
}

process RUN_FASTQC {
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
    cp ${adapter_file} \$HOME/.fastqc/adapters.txt
    
    # Run FastQC with adapter detection
    fastqc -q -t ${task.cpus} \\
        --adapters \$HOME/.fastqc/adapters.txt \\
        --extract \\
        --nogroup \\
        ${reads}
    """
}
