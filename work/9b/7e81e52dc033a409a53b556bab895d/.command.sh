#!/bin/bash -ue
# Create directory for adapter file
mkdir -p $HOME/.fastqc/
cp illumina_adapters.fa $HOME/.fastqc/adapters.fa

# Run FastQC with adapter detection
fastqc -q -t 2 \
    --adapters $HOME/.fastqc/adapters.fa \
    --extract \
    --nogroup \
    sample1_R1.fastq.gz sample1_R2.fastq.gz
