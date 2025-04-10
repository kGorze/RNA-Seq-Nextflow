#!/usr/bin/env nextflow

/*
 * Enhanced MultiQC module for RNA-Seq pipeline
 * Includes custom configuration and comprehensive reporting
 */

workflow MULTIQC {
    take:
    multiqc_files  // Channel with all QC files to be aggregated
    
    main:
    // Create MultiQC config
    CREATE_MULTIQC_CONFIG()
    config_file = CREATE_MULTIQC_CONFIG.out.config
    
    // Run MultiQC
    RUN_MULTIQC(multiqc_files, config_file)
    
    emit:
    report = RUN_MULTIQC.out.report
    data = RUN_MULTIQC.out.data
}

process CREATE_MULTIQC_CONFIG {
    label 'process_low'
    
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    output:
    path "multiqc_config.yaml", emit: config
    
    script:
    """
    cat > multiqc_config.yaml << EOL
# MultiQC Configuration File for RNA-Seq Pipeline

title: "RNA-Seq Analysis Report"
subtitle: "Quality Control and Processing Summary"
intro_text: "This report summarizes the quality control metrics and processing results for RNA-Seq data."

report_comment: "Generated by the Nextflow RNA-Seq Pipeline"

# Module-specific settings
fastqc_config:
    theoretical_gc: null

# Custom plot configuration
custom_plot_config:
    fastqc_per_base_sequence_quality:
        title: "Per Base Sequence Quality"
        xlab: "Position (bp)"
        ylab: "Quality Score"
    
    fastqc_per_sequence_quality_scores:
        title: "Per Sequence Quality Scores"
        xlab: "Quality Score"
        ylab: "Count"

# Output settings
output_dir: "${params.outdir}/multiqc"
make_data_dir: true
data_format: "tsv"

# Report appearance
template: "default"
EOL
    """
}

process RUN_MULTIQC {
    label 'process_low'
    
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path multiqc_files
    path config_file
    
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
    
    script:
    """
    # Run MultiQC with custom config
    multiqc -f \\
        -c ${config_file} \\
        -n multiqc_report.html \\
        -o . \\
        .
    """
}
