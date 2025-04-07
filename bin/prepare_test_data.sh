#!/bin/bash

# Script to download and prepare test dataset for RNA-Seq pipeline
# This script downloads a small subset of real RNA-Seq data from SRA
# and prepares it for testing the pipeline

# Create test data directory
mkdir -p test_data/fastq
cd test_data

# Download a small subset of RNA-Seq data from SRA
# Using E. coli RNA-Seq data (SRP013747) as it's small and well-characterized
echo "Downloading test data from SRA..."

# Download SRA toolkit if not already installed
if ! command -v fastq-dump &> /dev/null; then
    echo "Installing SRA toolkit..."
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf sratoolkit.current-ubuntu64.tar.gz
    export PATH=$PATH:$PWD/sratoolkit.*/bin
    rm sratoolkit.current-ubuntu64.tar.gz
fi

# Download a small subset of RNA-Seq data (2 samples, 2 conditions)
# Control samples
fastq-dump --split-files --gzip SRR515697 -O fastq
fastq-dump --split-files --gzip SRR515698 -O fastq

# Treatment samples
fastq-dump --split-files --gzip SRR515708 -O fastq
fastq-dump --split-files --gzip SRR515709 -O fastq

# Rename files to more intuitive names
mv fastq/SRR515697_1.fastq.gz fastq/control_rep1_R1.fastq.gz
mv fastq/SRR515697_2.fastq.gz fastq/control_rep1_R2.fastq.gz
mv fastq/SRR515698_1.fastq.gz fastq/control_rep2_R1.fastq.gz
mv fastq/SRR515698_2.fastq.gz fastq/control_rep2_R2.fastq.gz
mv fastq/SRR515708_1.fastq.gz fastq/treatment_rep1_R1.fastq.gz
mv fastq/SRR515708_2.fastq.gz fastq/treatment_rep1_R2.fastq.gz
mv fastq/SRR515709_1.fastq.gz fastq/treatment_rep2_R1.fastq.gz
mv fastq/SRR515709_2.fastq.gz fastq/treatment_rep2_R2.fastq.gz

# Download E. coli reference genome and annotation
echo "Downloading reference genome and annotation..."
wget -O genome.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-51/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz
wget -O annotation.gtf.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-51/gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.51.gtf.gz

# Decompress reference files
gunzip genome.fa.gz
gunzip annotation.gtf.gz

# Create design and contrasts files
echo "Creating experimental design files..."
cat > design.csv << EOL
sample_id,condition
control_rep1,control
control_rep2,control
treatment_rep1,treatment
treatment_rep2,treatment
EOL

cat > contrasts.csv << EOL
name,control,treatment
treatment_vs_control,control,treatment
EOL

# Create a README file for the test dataset
cat > README.md << EOL
# RNA-Seq Pipeline Test Dataset

This directory contains a small test dataset for the RNA-Seq Nextflow pipeline.

## Dataset Description

The test dataset consists of a subset of RNA-Seq data from *Escherichia coli* K-12 MG1655 
(SRA project SRP013747). The experiment compares gene expression between control and 
treatment conditions.

## Files

- **fastq/**: Contains paired-end FASTQ files for 4 samples (2 control, 2 treatment)
- **genome.fa**: *E. coli* K-12 MG1655 reference genome
- **annotation.gtf**: Gene annotation in GTF format
- **design.csv**: Sample metadata and experimental design
- **contrasts.csv**: Contrasts for differential expression analysis

## Usage

To run the pipeline with this test dataset:

\`\`\`bash
nextflow run main.nf \\
  -profile test,docker \\
  --reads 'test_data/fastq/*_R{1,2}.fastq.gz' \\
  --genome 'test_data/genome.fa' \\
  --gtf 'test_data/annotation.gtf' \\
  --design 'test_data/design.csv' \\
  --contrasts 'test_data/contrasts.csv' \\
  --outdir './results'
\`\`\`

## Source

The data was obtained from the Sequence Read Archive (SRA) project SRP013747:
https://www.ncbi.nlm.nih.gov/sra/?term=SRP013747
EOL

echo "Test dataset preparation complete!"
echo "The test dataset is available in the test_data directory."
