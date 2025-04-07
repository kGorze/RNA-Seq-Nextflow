#!/usr/bin/env nextflow

/*
 * Integrated reporting workflow for RNA-Seq pipeline
 * Combines all reporting modules into a cohesive reporting system
 */

// Import reporting modules
include { FINAL_REPORT_WORKFLOW } from './main'
include { PIPELINE_INFO_WORKFLOW } from './pipeline_info'
include { OUTPUT_ORGANIZATION_WORKFLOW } from './output_organization'

workflow REPORTING_WORKFLOW {
    take:
    fastqc_results    // Channel with FastQC results
    multiqc_results   // Channel with MultiQC results
    trimming_results  // Channel with trimming results
    alignment_results // Channel with alignment results
    alignment_stats   // Channel with alignment statistics
    quant_results     // Channel with quantification results
    quant_stats       // Channel with quantification statistics
    de_results        // Channel with differential expression results
    de_report         // Channel with differential expression report
    
    main:
    // Generate pipeline execution information
    PIPELINE_INFO_WORKFLOW()
    
    // Generate final report
    FINAL_REPORT_WORKFLOW(
        fastqc_results,
        multiqc_results.collect().ifEmpty([]),
        alignment_stats,
        quant_stats,
        de_results,
        de_report,
        PIPELINE_INFO_WORKFLOW.out.info
    )
    
    // Organize output files
    OUTPUT_ORGANIZATION_WORKFLOW(
        fastqc_results,
        multiqc_results,
        trimming_results,
        alignment_results,
        quant_results,
        de_results,
        FINAL_REPORT_WORKFLOW.out.report
    )
    
    emit:
    report = FINAL_REPORT_WORKFLOW.out.report
    summary = FINAL_REPORT_WORKFLOW.out.summary
    archive = FINAL_REPORT_WORKFLOW.out.archive
    documentation = OUTPUT_ORGANIZATION_WORKFLOW.out.documentation
    pipeline_info = PIPELINE_INFO_WORKFLOW.out.info
}

// Helper workflow to generate mock data for testing reporting
workflow MOCK_REPORTING_DATA {
    main:
    // Generate mock data for testing reporting
    GENERATE_MOCK_DATA()
    
    emit:
    fastqc = GENERATE_MOCK_DATA.out.fastqc
    multiqc = GENERATE_MOCK_DATA.out.multiqc
    alignment = GENERATE_MOCK_DATA.out.alignment
    alignment_stats = GENERATE_MOCK_DATA.out.alignment_stats
    quant = GENERATE_MOCK_DATA.out.quant
    quant_stats = GENERATE_MOCK_DATA.out.quant_stats
    de = GENERATE_MOCK_DATA.out.de
    de_report = GENERATE_MOCK_DATA.out.de_report
}

process GENERATE_MOCK_DATA {
    label 'process_low'
    
    publishDir "${params.outdir}/mock_data", mode: 'copy'
    
    output:
    path "fastqc/*", emit: fastqc
    path "multiqc/*", emit: multiqc
    path "alignment/*", emit: alignment
    path "alignment_stats.csv", emit: alignment_stats
    path "quant/*", emit: quant
    path "quant_stats.csv", emit: quant_stats
    path "de/*", emit: de
    path "de_report.html", emit: de_report
    
    script:
    """
    # Create directories
    mkdir -p fastqc multiqc alignment quant de
    
    # Create mock FastQC results
    echo "<html><body><h1>Mock FastQC Report</h1></body></html>" > fastqc/mock_fastqc.html
    
    # Create mock MultiQC report
    echo "<html><body><h1>Mock MultiQC Report</h1></body></html>" > multiqc/multiqc_report.html
    
    # Create mock alignment stats
    cat > alignment_stats.csv << EOL
Sample,Total_Reads,Mapped_Reads,Uniquely_Mapped,Multi_Mapped,Alignment_Rate
sample1,1000000,950000,900000,50000,95.0
sample2,1200000,1140000,1080000,60000,95.0
sample3,1100000,1045000,990000,55000,95.0
sample4,1300000,1235000,1170000,65000,95.0
EOL
    
    # Create mock quantification stats
    cat > quant_stats.csv << EOL
Sample,Total_Reads,Assigned_Reads,Detected_Genes,Assignment_Rate
sample1,950000,855000,15000,90.0
sample2,1140000,1026000,15200,90.0
sample3,1045000,940500,15100,90.0
sample4,1235000,1111500,15300,90.0
EOL
    
    # Create mock DE results
    cat > de/de_summary.csv << EOL
Contrast,Total_Genes,Significant_Genes,Upregulated,Downregulated
treatment_vs_control,20000,1500,800,700
EOL
    
    # Create mock DE report
    echo "<html><body><h1>Mock DE Report</h1></body></html>" > de_report.html
    
    # Create mock alignment files
    touch alignment/sample1.bam
    touch alignment/sample2.bam
    touch alignment/sample3.bam
    touch alignment/sample4.bam
    
    # Create mock quantification files
    touch quant/sample1.counts
    touch quant/sample2.counts
    touch quant/sample3.counts
    touch quant/sample4.counts
    
    # Create mock DE files
    touch de/treatment_vs_control_results.csv
    """
}
