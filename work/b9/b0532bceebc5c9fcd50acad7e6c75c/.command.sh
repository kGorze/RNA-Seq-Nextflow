#!/bin/bash -ue
trimmomatic PE -threads 2 \
    sample1_R1.fastq.gz sample1_R2.fastq.gz \
    sample1_1_trimmed.fastq.gz sample1_1_unpaired.fastq.gz \
    sample1_2_trimmed.fastq.gz sample1_2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 MINLEN:36 \
    2> sample1_trimming_report.txt
