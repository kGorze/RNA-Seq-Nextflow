#!/usr/bin/env nextflow

/*
 * Enhanced STAR alignment workflow for RNA-Seq pipeline
 * Includes genome indexing, alignment, and BAM processing
 */

workflow STAR {
    take:
        reads
        genome
        annotation

    main:
        STAR_INDEX(genome, annotation)
        STAR_ALIGN(reads, STAR_INDEX.out.index, annotation)

    emit:
        bam = STAR_ALIGN.out.bam
        counts = STAR_ALIGN.out.counts
        junctions = STAR_ALIGN.out.junctions
}

process STAR_INDEX {
    tag "${genome.simpleName}"
    container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    publishDir path: "${params.outdir}/star/index", mode: 'copy'
    
    memory '4 GB'
    time '1h'

    input:
        path genome
        path gtf

    output:
        path "star_index", emit: index

    script:
        """
        # Debug information
        echo "Running STAR indexing for ${genome.simpleName}"
        echo "Working directory: \$(pwd)"
        echo "Genome file: ${genome}"
        echo "GTF file: ${gtf}"
        
        # Create index directory
        mkdir -p star_index
        
        # Run STAR indexing
        STAR --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --genomeDir star_index \\
            --genomeFastaFiles ${genome} \\
            --sjdbGTFfile ${gtf} \\
            --sjdbOverhang 100 \\
            --limitGenomeGenerateRAM 4000000000
        
        # Verify the index was created
        echo "Contents of star_index directory:"
        ls -la star_index
        """
}

process STAR_ALIGN {
    tag "${meta.id}"
    container "quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    publishDir path: "${params.outdir}/star/aligned", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path index
        path gtf

    output:
        tuple val(meta), path("${meta.id}.Aligned.out.bam"), emit: bam
        tuple val(meta), path("${meta.id}.ReadsPerGene.out.tab"), emit: counts
        tuple val(meta), path("${meta.id}.SJ.out.tab"), emit: junctions

    script:
        def prefix = meta.id
        def read_args = reads instanceof List ? "--readFilesIn ${reads[0]} ${reads[1]}" : "--readFilesIn ${reads}"
        """
        # Debug information
        echo "Running STAR alignment for ${prefix}"
        echo "Working directory: \$(pwd)"
        echo "Contents of current directory:"
        ls -la
        
        # Create a truly temporary directory in the Docker container
        TMPDIR=\$(mktemp -d)
        echo "Created temporary directory: \$TMPDIR"
        chmod 777 \$TMPDIR
        
        # Run STAR with the container-local temporary directory
        STAR --runMode alignReads \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${index} \\
            --readFilesCommand zcat \\
            ${read_args} \\
            --outFileNamePrefix ${prefix}. \\
            --outTmpDir \$TMPDIR/star_tmp \\
            --outSAMtype BAM Unsorted \\
            --outSAMattributes Standard \\
            --outSAMunmapped Within \\
            --outFilterType BySJout \\
            --outFilterMultimapNmax 20 \\
            --outFilterMismatchNmax 999 \\
            --outFilterMismatchNoverReadLmax 0.04 \\
            --alignIntronMin 20 \\
            --alignIntronMax 1000000 \\
            --alignMatesGapMax 1000000 \\
            --alignSJoverhangMin 8 \\
            --alignSJDBoverhangMin 1 \\
            --sjdbGTFfile ${gtf} \\
            --quantMode GeneCounts
        """
}

process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/samtools:1.13--h8c37831_0'
    
    publishDir "${params.outdir}/star/sorted", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam
    path "${sample_id}.sorted.bam.bai", emit: bam_index
    
    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}

process SAMTOOLS_STATS {
    tag "$sample_id"
    label 'process_low'
    
    container 'quay.io/biocontainers/samtools:1.13--h8c37831_0'
    
    publishDir "${params.outdir}/star/stats", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path "${sample_id}.stats.txt", emit: stats
    path "${sample_id}.flagstat.txt", emit: flagstat
    path "${sample_id}.idxstats.txt", emit: idxstats
    
    script:
    """
    samtools stats ${bam} > ${sample_id}.stats.txt
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt
    """
}
