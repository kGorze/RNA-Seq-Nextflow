# RNA-Seq Pipeline Test Data

This directory contains test data generation scripts and sample data for the RNA-Seq pipeline.

## Test Data Generator

The `test_data_generator.py` script generates realistic test data for RNA-Seq analysis, including:

1. Paired-end FASTQ files with:
   - 10,000 reads per file
   - 100bp read length
   - Realistic quality scores following a normal distribution
   - Proper GC content distribution

2. Genome file (`test_genome.fa`):
   - 3 chromosomes
   - 100kb per chromosome
   - Realistic GC content distribution

3. Gene annotation file (`test_annotation.gtf`):
   - 50 protein-coding genes
   - Realistic gene structure including:
     - 5' and 3' UTRs
     - Multiple exons (2-5 per gene)
     - Introns
     - Proper strand information

## Generated Files

- `sample1_R1.fastq.gz`: Forward reads
- `sample1_R2.fastq.gz`: Reverse reads
- `test_genome.fa`: Reference genome
- `test_annotation.gtf`: Gene annotations

## Usage

To generate test data:

```bash
python test_data_generator.py
```

## Dependencies

- Python 3.x
- NumPy
- Biopython

## Data Structure

### FASTQ Files
- Paired-end reads with realistic quality scores
- Read length: 100bp
- Number of reads: 10,000 per file

### Genome File
- Format: FASTA
- Number of chromosomes: 3
- Chromosome length: 100kb each
- GC content: 40-60%

### Annotation File
- Format: GTF
- Features:
  - Genes
  - Exons
  - 5' UTRs
  - 3' UTRs
  - Introns
- Gene structure:
  - 2-5 exons per gene
  - 50-200bp UTRs
  - 50-500bp introns
  - 100-300bp exons 