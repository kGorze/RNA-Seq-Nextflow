#!/bin/bash

# Script to download a subset of real RNA-Seq data from public repositories
# This script downloads a small subset of human RNA-Seq data from the
# Sequence Read Archive (SRA) for testing the pipeline

# Create test data directory
mkdir -p test_data_human/fastq
cd test_data_human

# Download a small subset of human RNA-Seq data
# Using a subset of RNA-Seq data from human cell lines (GSE129240)
echo "Downloading human test data from SRA..."

# Download SRA toolkit if not already installed
if ! command -v fastq-dump &> /dev/null; then
    echo "Installing SRA toolkit..."
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf sratoolkit.current-ubuntu64.tar.gz
    export PATH=$PATH:$PWD/sratoolkit.*/bin
    rm sratoolkit.current-ubuntu64.tar.gz
fi

# Download a small subset of human RNA-Seq data (2 samples, 2 conditions)
# Control samples (HEK293 cells)
fastq-dump --split-files --gzip SRR8670432 -O fastq
fastq-dump --split-files --gzip SRR8670433 -O fastq

# Treatment samples (HEK293 cells with treatment)
fastq-dump --split-files --gzip SRR8670438 -O fastq
fastq-dump --split-files --gzip SRR8670439 -O fastq

# Rename files to more intuitive names
mv fastq/SRR8670432_1.fastq.gz fastq/control_rep1_R1.fastq.gz
mv fastq/SRR8670432_2.fastq.gz fastq/control_rep1_R2.fastq.gz
mv fastq/SRR8670433_1.fastq.gz fastq/control_rep2_R1.fastq.gz
mv fastq/SRR8670433_2.fastq.gz fastq/control_rep2_R2.fastq.gz
mv fastq/SRR8670438_1.fastq.gz fastq/treatment_rep1_R1.fastq.gz
mv fastq/SRR8670438_2.fastq.gz fastq/treatment_rep1_R2.fastq.gz
mv fastq/SRR8670439_1.fastq.gz fastq/treatment_rep2_R1.fastq.gz
mv fastq/SRR8670439_2.fastq.gz fastq/treatment_rep2_R2.fastq.gz

# Download a subset of human reference genome (chromosome 22) and annotation
echo "Downloading human reference genome (chr22) and annotation..."
wget -O genome.fa.gz http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
wget -O annotation.gtf.gz http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chromosome.22.gtf.gz

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
# Human RNA-Seq Pipeline Test Dataset

This directory contains a small test dataset of human RNA-Seq data for the RNA-Seq Nextflow pipeline.

## Dataset Description

The test dataset consists of a subset of RNA-Seq data from human HEK293 cells 
(GEO series GSE129240). The experiment compares gene expression between control and 
treatment conditions. For simplicity and to reduce computational requirements, only 
chromosome 22 is included in the reference genome.

## Files

- **fastq/**: Contains paired-end FASTQ files for 4 samples (2 control, 2 treatment)
- **genome.fa**: Human reference genome (chromosome 22 only)
- **annotation.gtf**: Gene annotation for chromosome 22 in GTF format
- **design.csv**: Sample metadata and experimental design
- **contrasts.csv**: Contrasts for differential expression analysis

## Usage

To run the pipeline with this test dataset:

\`\`\`bash
nextflow run main.nf \\
  -profile test,docker \\
  --reads 'test_data_human/fastq/*_R{1,2}.fastq.gz' \\
  --genome 'test_data_human/genome.fa' \\
  --gtf 'test_data_human/annotation.gtf' \\
  --design 'test_data_human/design.csv' \\
  --contrasts 'test_data_human/contrasts.csv' \\
  --outdir './results'
\`\`\`

## Source

The data was obtained from the Gene Expression Omnibus (GEO) series GSE129240:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129240
EOL

echo "Human test dataset preparation complete!"
echo "The human test dataset is available in the test_data_human directory."
