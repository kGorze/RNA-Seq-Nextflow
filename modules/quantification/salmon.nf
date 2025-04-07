#!/usr/bin/env nextflow

/*
 * Enhanced Salmon quantification workflow for RNA-Seq pipeline
 * Includes transcriptome generation, indexing, and quantification
 */

workflow SALMON_WORKFLOW {
    take:
    reads_ch     // Channel with sample_id and reads
    genome_file  // Reference genome FASTA
    gtf_file     // Gene annotation GTF
    
    main:
    // Generate transcriptome from genome and annotation
    GENERATE_TRANSCRIPTOME(genome_file, gtf_file)
    
    // Build Salmon index
    SALMON_INDEX(GENERATE_TRANSCRIPTOME.out.transcriptome)
    
    // Quantify with Salmon
    SALMON_QUANT(reads_ch, SALMON_INDEX.out.index)
    
    // Generate gene-level counts from transcript-level counts
    SALMON_TXIMPORT(SALMON_QUANT.out.results.collect(), gtf_file)
    
    emit:
    counts = SALMON_QUANT.out.results
    merged_counts = SALMON_TXIMPORT.out.gene_counts
    tx_counts = SALMON_TXIMPORT.out.tx_counts
    logs = SALMON_QUANT.out.logs
}

process GENERATE_TRANSCRIPTOME {
    label 'process_medium'
    
    container 'quay.io/biocontainers/gffread:0.12.7--h9ee0642_0'
    
    publishDir "${params.outdir}/salmon/transcriptome", mode: 'copy'
    
    input:
    path genome
    path gtf
    
    output:
    path "transcripts.fa", emit: transcriptome
    
    script:
    """
    gffread ${gtf} -g ${genome} -w transcripts.fa
    """
}

process SALMON_INDEX {
    label 'process_high'
    
    container 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
    
    publishDir "${params.outdir}/salmon/index", mode: 'copy'
    
    input:
    path transcriptome
    
    output:
    path "salmon_index", emit: index
    
    script:
    """
    salmon index \\
        -p ${task.cpus} \\
        -t ${transcriptome} \\
        -i salmon_index \\
        --gencode
    """
}

process SALMON_QUANT {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
    
    publishDir "${params.outdir}/salmon/quant", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path index
    
    output:
    tuple val(sample_id), path("${sample_id}"), emit: results
    path "${sample_id}/logs/*", emit: logs
    
    script:
    def prefix = "${sample_id}"
    
    // Determine library type based on strandedness parameter
    def lib_type = ""
    if (params.paired_end) {
        lib_type = params.strandedness == 'forward' ? 'ISF' : 
                   params.strandedness == 'reverse' ? 'ISR' : 'IU'
    } else {
        lib_type = params.strandedness == 'forward' ? 'SF' : 
                   params.strandedness == 'reverse' ? 'SR' : 'U'
    }
    
    // Determine read arguments based on paired-end parameter
    def read_args = ""
    if (params.paired_end) {
        read_args = "-1 ${reads[0]} -2 ${reads[1]}"
    } else {
        read_args = "-r ${reads}"
    }
    
    """
    salmon quant \\
        -p ${task.cpus} \\
        -i ${index} \\
        -l ${lib_type} \\
        ${read_args} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        -o ${prefix}
    """
}

process SALMON_TXIMPORT {
    label 'process_medium'
    
    container 'quay.io/biocontainers/bioconductor-tximport:1.20.0--r41hdfd78af_0'
    
    publishDir "${params.outdir}/salmon/tximport", mode: 'copy'
    
    input:
    path quant_dirs
    path gtf
    
    output:
    path "gene_counts.csv", emit: gene_counts
    path "transcript_counts.csv", emit: tx_counts
    path "tximport_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    library(tximport)
    library(GenomicFeatures)
    
    # Get list of quant.sf files
    quant_files <- list.files(path=".", pattern="quant.sf", recursive=TRUE, full.names=TRUE)
    names(quant_files) <- basename(dirname(quant_files))
    
    # Create tx2gene mapping from GTF
    txdb <- GenomicFeatures::makeTxDbFromGFF("${gtf}")
    k <- keys(txdb, keytype="TXNAME")
    tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
    
    # Import transcript-level estimates
    txi <- tximport(quant_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
    
    # Write gene-level counts to file
    write.csv(txi\$counts, file="gene_counts.csv")
    
    # Write transcript-level counts to file
    txi_tx <- tximport(quant_files, type="salmon", txOut=TRUE)
    write.csv(txi_tx\$counts, file="transcript_counts.csv")
    
    # Write summary
    sink("tximport_summary.txt")
    cat("tximport summary\\n")
    cat("----------------\\n")
    cat("Number of samples:", length(quant_files), "\\n")
    cat("Number of genes:", nrow(txi\$counts), "\\n")
    cat("Number of transcripts:", nrow(txi_tx\$counts), "\\n")
    sink()
    """
}
