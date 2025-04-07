#!/bin/bash -ue
mkdir -p star_index

STAR --runMode genomeGenerate \
    --runThreadN 2 \
    --genomeDir star_index \
    --genomeFastaFiles test_genome.fa \
    --sjdbGTFfile test_annotation.gtf \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11
