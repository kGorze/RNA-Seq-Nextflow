#!/bin/bash

# Script to download and subsample public RNA-Seq data for testing
# This script creates a very small subset of RNA-Seq data for quick testing

# Create test data directory
mkdir -p test_data_mini/fastq
cd test_data_mini

# Download a small subset of RNA-Seq data and subsample it further
echo "Downloading and subsampling RNA-Seq data..."

# Download SRA toolkit if not already installed
if ! command -v fastq-dump &> /dev/null; then
    echo "Installing SRA toolkit..."
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf sratoolkit.current-ubuntu64.tar.gz
    export PATH=$PATH:$PWD/sratoolkit.*/bin
    rm sratoolkit.current-ubuntu64.tar.gz
fi

# Download a single sample from SRA
fastq-dump --split-files --gzip SRR515697 -O fastq_temp

# Install seqtk for subsampling if not already installed
if ! command -v seqtk &> /dev/null; then
    echo "Installing seqtk..."
    git clone https://github.com/lh3/seqtk.git
    cd seqtk
    make
    export PATH=$PATH:$PWD
    cd ..
fi

# Subsample to create very small FASTQ files (1000 reads per sample)
mkdir -p fastq
seqtk sample -s100 fastq_temp/SRR515697_1.fastq.gz 1000 | gzip > fastq/control_rep1_R1.fastq.gz
seqtk sample -s100 fastq_temp/SRR515697_2.fastq.gz 1000 | gzip > fastq/control_rep1_R2.fastq.gz
seqtk sample -s200 fastq_temp/SRR515697_1.fastq.gz 1000 | gzip > fastq/control_rep2_R1.fastq.gz
seqtk sample -s200 fastq_temp/SRR515697_2.fastq.gz 1000 | gzip > fastq/control_rep2_R2.fastq.gz
seqtk sample -s300 fastq_temp/SRR515697_1.fastq.gz 1000 | gzip > fastq/treatment_rep1_R1.fastq.gz
seqtk sample -s300 fastq_temp/SRR515697_2.fastq.gz 1000 | gzip > fastq/treatment_rep1_R2.fastq.gz
seqtk sample -s400 fastq_temp/SRR515697_1.fastq.gz 1000 | gzip > fastq/treatment_rep2_R1.fastq.gz
seqtk sample -s400 fastq_temp/SRR515697_2.fastq.gz 1000 | gzip > fastq/treatment_rep2_R2.fastq.gz

# Clean up temporary files
rm -rf fastq_temp

# Create a very small reference genome and annotation
echo "Creating mini reference genome and annotation..."

# Create a mini reference genome with just a few genes
cat > genome.fa << EOL
>mini_chr1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC
>mini_chr2
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC
EOL

# Create a mini annotation
cat > annotation.gtf << EOL
mini_chr1	mini	gene	1	100	.	+	.	gene_id "gene1"; gene_name "Gene1";
mini_chr1	mini	transcript	1	100	.	+	.	gene_id "gene1"; transcript_id "tx1"; gene_name "Gene1";
mini_chr1	mini	exon	1	50	.	+	.	gene_id "gene1"; transcript_id "tx1"; exon_id "exon1"; gene_name "Gene1";
mini_chr1	mini	exon	51	100	.	+	.	gene_id "gene1"; transcript_id "tx1"; exon_id "exon2"; gene_name "Gene1";
mini_chr1	mini	gene	150	250	.	-	.	gene_id "gene2"; gene_name "Gene2";
mini_chr1	mini	transcript	150	250	.	-	.	gene_id "gene2"; transcript_id "tx2"; gene_name "Gene2";
mini_chr1	mini	exon	150	250	.	-	.	gene_id "gene2"; transcript_id "tx2"; exon_id "exon3"; gene_name "Gene2";
mini_chr2	mini	gene	1	100	.	+	.	gene_id "gene3"; gene_name "Gene3";
mini_chr2	mini	transcript	1	100	.	+	.	gene_id "gene3"; transcript_id "tx3"; gene_name "Gene3";
mini_chr2	mini	exon	1	100	.	+	.	gene_id "gene3"; transcript_id "tx3"; exon_id "exon4"; gene_name "Gene3";
mini_chr2	mini	gene	150	250	.	-	.	gene_id "gene4"; gene_name "Gene4";
mini_chr2	mini	transcript	150	250	.	-	.	gene_id "gene4"; transcript_id "tx4"; gene_name "Gene4";
mini_chr2	mini	exon	150	200	.	-	.	gene_id "gene4"; transcript_id "tx4"; exon_id "exon5"; gene_name "Gene4";
mini_chr2	mini	exon	201	250	.	-	.	gene_id "gene4"; transcript_id "tx4"; exon_id "exon6"; gene_name "Gene4";
EOL

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
# Mini RNA-Seq Pipeline Test Dataset

This directory contains a minimal test dataset for the RNA-Seq Nextflow pipeline.
This dataset is designed for very quick testing and CI/CD purposes.

## Dataset Description

The test dataset consists of a tiny subset of RNA-Seq data (1000 reads per sample)
and a minimal synthetic reference genome with just a few genes.

## Files

- **fastq/**: Contains paired-end FASTQ files for 4 samples (2 control, 2 treatment)
- **genome.fa**: Minimal synthetic reference genome
- **annotation.gtf**: Minimal gene annotation in GTF format
- **design.csv**: Sample metadata and experimental design
- **contrasts.csv**: Contrasts for differential expression analysis

## Usage

To run the pipeline with this test dataset:

\`\`\`bash
nextflow run main.nf \\
  -profile test,docker \\
  --reads 'test_data_mini/fastq/*_R{1,2}.fastq.gz' \\
  --genome 'test_data_mini/genome.fa' \\
  --gtf 'test_data_mini/annotation.gtf' \\
  --design 'test_data_mini/design.csv' \\
  --contrasts 'test_data_mini/contrasts.csv' \\
  --outdir './results'
\`\`\`

## Note

This is a minimal dataset for quick testing only. For more realistic testing,
use the other test datasets provided in the repository.
EOL

echo "Mini test dataset preparation complete!"
echo "The mini test dataset is available in the test_data_mini directory."
