#!/usr/bin/env nextflow

/*
 * RNA-Seq Nextflow Pipeline
 * A reproducible RNA-Seq analysis pipeline implemented in Nextflow
 */

// Enable DSL2
nextflow.enable.dsl = 2

// Import modules
include { FASTQC } from './modules/fastqc/main'
include { MULTIQC } from './modules/multiqc/main'
include { TRIMMOMATIC } from './modules/trimming/main'
include { STAR } from './modules/alignment/star'
include { HISAT2 } from './modules/alignment/hisat2'
include { SALMON } from './modules/quantification/salmon'
include { FEATURECOUNTS } from './modules/quantification/featurecounts'
include { DESEQ2 } from './modules/deseq2/main'

// Print pipeline header
log.info """
=======================================================
RNA-Seq Nextflow Pipeline v${workflow.manifest.version}
=======================================================
Reads        : ${params.reads}
Genome       : ${params.genome}
Annotation   : ${params.gtf}
Aligner      : ${params.aligner}
Output dir   : ${params.outdir}
=======================================================
"""

// Check mandatory parameters
if (params.reads == null) {
    exit 1, "Input reads not specified. Please provide --reads parameter."
}
if (params.genome == null) {
    exit 1, "Reference genome not specified. Please provide --genome parameter."
}
if (params.gtf == null) {
    exit 1, "Gene annotation not specified. Please provide --gtf parameter."
}

// Define workflow
workflow {
    // Channel for input reads
    if (params.paired_end) {
        Channel
            .fromFilePairs(params.reads, checkIfExists: true)
            .set { read_pairs_ch }
    } else {
        Channel
            .fromPath(params.reads, checkIfExists: true)
            .map { file -> tuple(file.simpleName, file) }
            .set { read_pairs_ch }
    }
    
    // Reference genome and annotation
    genome_file = file(params.genome)
    gtf_file = file(params.gtf)
    
    // Quality control
    if (!params.skip_qc) {
        FASTQC(read_pairs_ch)
    }
    
    // Trimming
    if (!params.skip_trimming) {
        TRIMMOMATIC(read_pairs_ch)
        trimmed_reads_ch = TRIMMOMATIC.out.trimmed_reads
        trimming_logs = TRIMMOMATIC.out.log
    } else {
        trimmed_reads_ch = read_pairs_ch
    }
    
    // Alignment and quantification
    if (params.aligner == 'star') {
        STAR(trimmed_reads_ch, genome_file, gtf_file)
        FEATURECOUNTS(STAR.out.bam, gtf_file)
        counts_ch = FEATURECOUNTS.out.merged_gene_counts
    } else if (params.aligner == 'hisat2') {
        HISAT2(trimmed_reads_ch, genome_file, gtf_file)
        FEATURECOUNTS(HISAT2.out.bam, gtf_file)
        counts_ch = FEATURECOUNTS.out.merged_gene_counts
    } else if (params.aligner == 'salmon') {
        SALMON(trimmed_reads_ch, genome_file, gtf_file)
        counts_ch = SALMON.out.counts
    }
    
    // Differential expression analysis
    if (params.design != null && params.contrasts != null) {
        design_file = file(params.design)
        contrasts_file = file(params.contrasts)
        DESEQ2(counts_ch, design_file, contrasts_file)
    }
    
    // Quality control and MultiQC report
    if (!params.skip_qc) {
        // Collect all QC files using Channel.fromPath after processes have run
        Channel
            .fromPath("${params.outdir}/**{_fastqc.zip,_fastqc.html,Log.*,*_stats.txt,*.summary}", hidden: true)
            .collect()
            .set { multiqc_files }
        
        // Generate MultiQC report
        MULTIQC(multiqc_files)
    }
}

// Workflow completion notification
workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'Failed'}"
    log.info "Execution duration: ${workflow.duration}"
}
