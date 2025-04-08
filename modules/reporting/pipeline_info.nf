#!/usr/bin/env nextflow

/*
 * Pipeline execution information module for RNA-Seq pipeline
 * Captures and reports pipeline execution metrics
 */

workflow PIPELINE_INFO_WORKFLOW {
    main:
    // Generate pipeline execution information
    GENERATE_PIPELINE_INFO()
    
    emit:
    info = GENERATE_PIPELINE_INFO.out.info
}

process GENERATE_PIPELINE_INFO {
    label 'process_low'
    
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    
    output:
    path "pipeline_info.txt", emit: info
    path "execution_trace.txt", emit: trace
    path "software_versions.yml", emit: versions
    
    script:
    """
    # Generate pipeline information
    cat > pipeline_info.txt << EOL
Pipeline version: ${workflow.manifest.version}
Run name: ${workflow.runName}
Run date: \$(date)
Execution time: ${workflow.duration}
Command line: ${workflow.commandLine}
Container: ${workflow.container}
Profile: ${workflow.profile}
EOL
    
    # Generate execution trace
    cat > execution_trace.txt << EOL
Workflow session: ${workflow.sessionId}
Workflow start: ${workflow.start}
Workflow complete: ${workflow.complete}
Workflow duration: ${workflow.duration}
Workflow success: ${workflow.success}
Workflow exit status: ${workflow.exitStatus}
EOL
    
    # Generate software versions
    cat > software_versions.yml << EOL
---
nextflow: ${nextflow.version}
fastqc: \$(fastqc --version 2>&1 | grep -o 'FastQC v[0-9.]*' | sed 's/FastQC v//') || 'not installed'
multiqc: \$(multiqc --version 2>&1 | grep -o 'multiqc, version [0-9.]*' | sed 's/multiqc, version //') || 'not installed'
star: \$(STAR --version 2>&1 | grep -o 'STAR_[0-9.]*' | sed 's/STAR_//') || 'not installed'
hisat2: \$(hisat2 --version 2>&1 | grep -o 'version [0-9.]*' | sed 's/version //') || 'not installed'
salmon: \$(salmon --version 2>&1 | grep -o 'salmon [0-9.]*' | sed 's/salmon //') || 'not installed'
featurecounts: \$(featureCounts -v 2>&1 | grep -o 'featureCounts v[0-9.]*' | sed 's/featureCounts v//') || 'not installed'
r: \$(R --version 2>&1 | grep -o 'R version [0-9.]*' | sed 's/R version //') || 'not installed'
EOL
    """
}
