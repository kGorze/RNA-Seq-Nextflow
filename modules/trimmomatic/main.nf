process CREATE_ADAPTER_FILE {
    container "ubuntu:20.04"
    
    output:
    path "TruSeq3-PE.fa", emit: adapter_file
    
    script:
    """
    cat > TruSeq3-PE.fa << EOL
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
EOL
    """
}

process RUN_TRIMMOMATIC {
    container "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    
    input:
    tuple val(sample), path(reads)
    path adapters
    
    output:
    tuple val(sample), path("${sample}_paired_R*.fastq.gz"), emit: trimmed_reads
    path "trimmomatic.log", emit: log_file
    
    script:
    """
    trimmomatic PE \\
        -threads 1 \\
        -phred33 \\
        ${reads[0]} \\
        ${reads[1]} \\
        ${sample}_paired_R1.fastq.gz \\
        ${sample}_unpaired_R1.fastq.gz \\
        ${sample}_paired_R2.fastq.gz \\
        ${sample}_unpaired_R2.fastq.gz \\
        ILLUMINACLIP:${adapters}:2:30:10 \\
        LEADING:3 \\
        TRAILING:3 \\
        SLIDINGWINDOW:4:15 \\
        MINLEN:36 \\
        2>&1 | tee trimmomatic.log
    """
}

workflow TRIMMOMATIC {
    take:
    reads
    
    main:
    adapter_file = CREATE_ADAPTER_FILE()
    RUN_TRIMMOMATIC(reads, adapter_file.adapter_file)
    
    emit:
    trimmed_reads = RUN_TRIMMOMATIC.out.trimmed_reads
    log_file = RUN_TRIMMOMATIC.out.log_file
} 