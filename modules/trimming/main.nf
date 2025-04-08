#!/usr/bin/env nextflow

/*
 * Trimmomatic module for RNA-Seq pipeline
 */

workflow TRIMMOMATIC {
    take:
    reads_ch  // Channel with sample_id and reads

    main:
    // Create adapter files
    CREATE_ADAPTER_FILES()
    
    // Run Trimmomatic
    RUN_TRIMMOMATIC(reads_ch, CREATE_ADAPTER_FILES.out.pe_adapters, CREATE_ADAPTER_FILES.out.se_adapters)
    
    emit:
    trimmed_reads = RUN_TRIMMOMATIC.out.trimmed_reads
    log = RUN_TRIMMOMATIC.out.log
}

process CREATE_ADAPTER_FILES {
    label 'process_low'
    
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    output:
    path "TruSeq3-PE.fa", emit: pe_adapters
    path "TruSeq3-SE.fa", emit: se_adapters
    
    script:
    """
    # Create PE adapter file
    cat > TruSeq3-PE.fa << EOL
>TruSeq_Universal_Adapter
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq_Adapter_Index_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_2
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_3
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_4
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_5
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_6
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_7
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_8
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
EOL

    # Create SE adapter file
    cat > TruSeq3-SE.fa << EOL
>TruSeq_Universal_Adapter
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq_Adapter_Index_1
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_2
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_3
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_4
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_5
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_6
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_7
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_Adapter_Index_8
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
EOL
    """
}

process RUN_TRIMMOMATIC {
    tag "$sample_id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
    
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path pe_adapters
    path se_adapters
    
    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: log
    
    script:
    def prefix = "${sample_id}"
    if (params.paired_end) {
        """
        trimmomatic PE -threads ${task.cpus} \\
            ${reads[0]} ${reads[1]} \\
            ${prefix}_1_trimmed.fastq.gz ${prefix}_1_unpaired.fastq.gz \\
            ${prefix}_2_trimmed.fastq.gz ${prefix}_2_unpaired.fastq.gz \\
            ILLUMINACLIP:${pe_adapters}:2:30:10:2:keepBothReads \\
            LEADING:3 TRAILING:3 MINLEN:36 \\
            2> ${prefix}_trimming_report.txt
        """
    } else {
        """
        trimmomatic SE -threads ${task.cpus} \\
            ${reads} \\
            ${prefix}_trimmed.fastq.gz \\
            ILLUMINACLIP:${se_adapters}:2:30:10 \\
            LEADING:3 TRAILING:3 MINLEN:36 \\
            2> ${prefix}_trimming_report.txt
        """
    }
}
