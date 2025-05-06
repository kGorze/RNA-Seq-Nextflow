#!/usr/bin/env nextflow

/*
 * Quality Control workflow for RNA-Seq pipeline
 * Integrates FastQC, MultiQC, and QC reporting
 */

// Import required modules
include { FASTQC_WORKFLOW } from './fastqc/main'
include { MULTIQC_WORKFLOW } from './multiqc/main'
include { QC_REPORT } from './fastqc/qc_report'
include { GENERATE_TEST_DATA } from './fastqc/test_data'

workflow QC_WORKFLOW {
    take:
    reads_ch  // Channel with sample_id and reads
    
    main:
    // Run FastQC workflow
    FASTQC_WORKFLOW(reads_ch)
    
    // Collect FastQC outputs for MultiQC
    fastqc_outputs = FASTQC_WORKFLOW.out.zip.collect()
    
    // Run MultiQC workflow
    MULTIQC_WORKFLOW(fastqc_outputs)
    
    // Generate QC report
    QC_REPORT(
        FASTQC_WORKFLOW.out.zip.collect(),
        MULTIQC_WORKFLOW.out.data
    )
    
    emit:
    fastqc_zip = FASTQC_WORKFLOW.out.zip
    fastqc_html = FASTQC_WORKFLOW.out.html
    multiqc_report = MULTIQC_WORKFLOW.out.report
    qc_report = QC_REPORT.out.report
}

// Workflow for generating test data
workflow GENERATE_TEST_DATASET {
    main:
    GENERATE_TEST_DATA()
    
    emit:
    reads = GENERATE_TEST_DATA.out.reads
    genome = GENERATE_TEST_DATA.out.genome
    gtf = GENERATE_TEST_DATA.out.gtf
}
