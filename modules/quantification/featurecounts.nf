#!/usr/bin/env nextflow

/*
 * Enhanced featureCounts quantification workflow for RNA-Seq pipeline
 * Includes gene-level and transcript-level counting
 */

workflow FEATURECOUNTS_WORKFLOW {
    take:
    bam_ch      // Channel with sample_id and BAM files
    gtf_file    // Gene annotation GTF
    
    main:
    // Gene-level counting
    FEATURECOUNTS_GENE(bam_ch, gtf_file)
    
    // Transcript-level counting
    FEATURECOUNTS_TRANSCRIPT(bam_ch, gtf_file)
    
    // Merge counts into a single matrix
    MERGE_COUNTS(
        FEATURECOUNTS_GENE.out.counts.collect(),
        FEATURECOUNTS_TRANSCRIPT.out.counts.collect()
    )
    
    emit:
    gene_counts = FEATURECOUNTS_GENE.out.counts
    tx_counts = FEATURECOUNTS_TRANSCRIPT.out.counts
    merged_gene_counts = MERGE_COUNTS.out.merged_gene_counts
    merged_tx_counts = MERGE_COUNTS.out.merged_tx_counts
    gene_logs = FEATURECOUNTS_GENE.out.log
    tx_logs = FEATURECOUNTS_TRANSCRIPT.out.log
}

process FEATURECOUNTS_GENE {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0'
    
    publishDir "${params.outdir}/featurecounts/gene", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path gtf
    
    output:
    tuple val(sample_id), path("${sample_id}.gene_counts.txt"), emit: counts
    path "${sample_id}.gene_counts.txt.summary", emit: log
    
    script:
    def prefix = "${sample_id}"
    
    // Set strand parameter based on strandedness
    def strand_param = params.strandedness == 'forward' ? '-s 1' : 
                       params.strandedness == 'reverse' ? '-s 2' : '-s 0'
    
    // Set paired-end parameter if needed
    def paired_end = params.paired_end ? '-p' : ''
    
    """
    featureCounts \\
        -T ${task.cpus} \\
        ${paired_end} \\
        ${strand_param} \\
        -a ${gtf} \\
        -g gene_id \\
        -o ${prefix}.gene_counts.txt \\
        ${bam}
    """
}

process FEATURECOUNTS_TRANSCRIPT {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0'
    
    publishDir "${params.outdir}/featurecounts/transcript", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path gtf
    
    output:
    tuple val(sample_id), path("${sample_id}.transcript_counts.txt"), emit: counts
    path "${sample_id}.transcript_counts.txt.summary", emit: log
    
    script:
    def prefix = "${sample_id}"
    
    // Set strand parameter based on strandedness
    def strand_param = params.strandedness == 'forward' ? '-s 1' : 
                       params.strandedness == 'reverse' ? '-s 2' : '-s 0'
    
    // Set paired-end parameter if needed
    def paired_end = params.paired_end ? '-p' : ''
    
    """
    featureCounts \\
        -T ${task.cpus} \\
        ${paired_end} \\
        ${strand_param} \\
        -a ${gtf} \\
        -g transcript_id \\
        -o ${prefix}.transcript_counts.txt \\
        ${bam}
    """
}

process MERGE_COUNTS {
    label 'process_low'
    
    container 'quay.io/biocontainers/r-base:4.1.0'
    
    publishDir "${params.outdir}/featurecounts/merged", mode: 'copy'
    
    input:
    path gene_counts
    path tx_counts
    
    output:
    path "merged_gene_counts.csv", emit: merged_gene_counts
    path "merged_transcript_counts.csv", emit: merged_tx_counts
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Function to merge count files
    merge_count_files <- function(count_files, output_file) {
        # Initialize empty data frame for merged counts
        merged_counts <- NULL
        sample_names <- c()
        
        # Process each count file
        for (file in count_files) {
            # Extract sample name from filename
            sample_name <- sub("\\..*counts\\.txt\$", "", basename(file))
            sample_names <- c(sample_names, sample_name)
            
            # Read count data (skip header lines)
            counts <- read.table(file, header=TRUE, skip=1)
            
            if (is.null(merged_counts)) {
                # First file - initialize with gene IDs and lengths
                merged_counts <- data.frame(
                    gene_id = counts[,1],
                    length = counts\$Length
                )
                rownames(merged_counts) <- counts[,1]
            }
            
            # Add counts from this sample
            merged_counts[[sample_name]] <- counts\$count
        }
        
        # Write merged counts to file
        write.csv(merged_counts, file=output_file, row.names=FALSE)
        
        return(sample_names)
    }
    
    # Get list of gene count files
    gene_files <- list.files(path=".", pattern="gene_counts.txt\$", full.names=TRUE)
    
    # Get list of transcript count files
    tx_files <- list.files(path=".", pattern="transcript_counts.txt\$", full.names=TRUE)
    
    # Merge gene counts
    gene_samples <- merge_count_files(gene_files, "merged_gene_counts.csv")
    
    # Merge transcript counts
    tx_samples <- merge_count_files(tx_files, "merged_transcript_counts.csv")
    
    # Print summary
    cat("Merged counts from", length(gene_samples), "samples\\n")
    cat("Samples:", paste(gene_samples, collapse=", "), "\\n")
    """
}
