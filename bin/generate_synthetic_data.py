#!/usr/bin/env python3

"""
Script to generate synthetic test data for RNA-Seq pipeline
This creates a small synthetic dataset that can be used for quick testing
"""

import os
import gzip
import random
import argparse
from pathlib import Path

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Generate synthetic RNA-Seq test data')
    parser.add_argument('--output_dir', type=str, default='test_data_synthetic',
                        help='Output directory for test data')
    parser.add_argument('--num_samples', type=int, default=4,
                        help='Number of samples to generate (half control, half treatment)')
    parser.add_argument('--num_reads', type=int, default=10000,
                        help='Number of reads per sample')
    parser.add_argument('--read_length', type=int, default=100,
                        help='Length of reads')
    parser.add_argument('--genome_size', type=int, default=50000,
                        help='Size of synthetic genome')
    parser.add_argument('--num_genes', type=int, default=200,
                        help='Number of genes in synthetic genome')
    return parser.parse_args()

def random_dna(length):
    """Generate random DNA sequence"""
    return ''.join(random.choice('ACGT') for _ in range(length))

def generate_fastq_entry(read_id, seq, qual_score=30):
    """Generate a FASTQ entry"""
    quality = chr(33 + qual_score) * len(seq)
    return f"@{read_id}\n{seq}\n+\n{quality}\n"

def generate_fastq_files(output_dir, num_samples, num_reads, read_length):
    """Generate paired-end FASTQ files"""
    fastq_dir = os.path.join(output_dir, 'fastq')
    os.makedirs(fastq_dir, exist_ok=True)
    
    # Generate half control, half treatment samples
    num_control = num_samples // 2
    num_treatment = num_samples - num_control
    
    # Generate control samples
    for i in range(1, num_control + 1):
        r1_file = os.path.join(fastq_dir, f"control_rep{i}_R1.fastq.gz")
        r2_file = os.path.join(fastq_dir, f"control_rep{i}_R2.fastq.gz")
        
        with gzip.open(r1_file, 'wt') as r1, gzip.open(r2_file, 'wt') as r2:
            for j in range(num_reads):
                read_id = f"control_rep{i}_read{j}"
                seq1 = random_dna(read_length)
                seq2 = random_dna(read_length)
                r1.write(generate_fastq_entry(read_id, seq1))
                r2.write(generate_fastq_entry(read_id, seq2))
    
    # Generate treatment samples
    for i in range(1, num_treatment + 1):
        r1_file = os.path.join(fastq_dir, f"treatment_rep{i}_R1.fastq.gz")
        r2_file = os.path.join(fastq_dir, f"treatment_rep{i}_R2.fastq.gz")
        
        with gzip.open(r1_file, 'wt') as r1, gzip.open(r2_file, 'wt') as r2:
            for j in range(num_reads):
                read_id = f"treatment_rep{i}_read{j}"
                seq1 = random_dna(read_length)
                seq2 = random_dna(read_length)
                r1.write(generate_fastq_entry(read_id, seq1))
                r2.write(generate_fastq_entry(read_id, seq2))
    
    return [f"control_rep{i}" for i in range(1, num_control + 1)] + \
           [f"treatment_rep{i}" for i in range(1, num_treatment + 1)]

def generate_genome(output_dir, genome_size):
    """Generate a synthetic reference genome"""
    genome_file = os.path.join(output_dir, 'genome.fa')
    
    with open(genome_file, 'w') as f:
        # Create a single chromosome
        f.write(">chr1\n")
        genome_seq = random_dna(genome_size)
        
        # Write sequence in lines of 80 characters
        for i in range(0, len(genome_seq), 80):
            f.write(genome_seq[i:i+80] + "\n")
    
    return genome_seq

def generate_annotation(output_dir, genome_seq, num_genes):
    """Generate a synthetic gene annotation"""
    gtf_file = os.path.join(output_dir, 'annotation.gtf')
    
    with open(gtf_file, 'w') as f:
        # Write header
        f.write('##description: Synthetic annotation for RNA-Seq pipeline testing\n')
        f.write('#!genome-build Synthetic\n')
        
        genome_length = len(genome_seq)
        gene_length = genome_length // num_genes
        
        # Create genes with transcripts and exons
        for gene_id in range(1, num_genes + 1):
            # Determine gene boundaries
            gene_start = (gene_id - 1) * gene_length + 1
            gene_end = gene_start + gene_length - 1
            if gene_end > genome_length:
                gene_end = genome_length
            
            strand = '+' if random.random() > 0.5 else '-'
            
            # Gene entry
            f.write(f'chr1\tsynthetic\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t')
            f.write(f'gene_id "gene{gene_id}"; gene_name "Gene{gene_id}";\n')
            
            # Create 1-2 transcripts for each gene
            num_transcripts = random.randint(1, 2)
            for tx_id in range(1, num_transcripts + 1):
                transcript_id = f"transcript{gene_id}_{tx_id}"
                tx_start = gene_start
                tx_end = gene_end
                
                # Transcript entry
                f.write(f'chr1\tsynthetic\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t')
                f.write(f'gene_id "gene{gene_id}"; transcript_id "{transcript_id}"; gene_name "Gene{gene_id}";\n')
                
                # Create 2-5 exons for each transcript
                num_exons = random.randint(2, 5)
                exon_length = (tx_end - tx_start + 1) // num_exons
                
                for exon_id in range(1, num_exons + 1):
                    exon_start = tx_start + (exon_id - 1) * exon_length
                    exon_end = exon_start + exon_length - 1
                    
                    if exon_id == num_exons:
                        exon_end = tx_end
                    
                    # Exon entry
                    f.write(f'chr1\tsynthetic\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t')
                    f.write(f'gene_id "gene{gene_id}"; transcript_id "{transcript_id}"; exon_number "{exon_id}"; gene_name "Gene{gene_id}";\n')

def generate_design_files(output_dir, sample_ids):
    """Generate design and contrasts files"""
    design_file = os.path.join(output_dir, 'design.csv')
    contrasts_file = os.path.join(output_dir, 'contrasts.csv')
    
    # Create design file
    with open(design_file, 'w') as f:
        f.write("sample_id,condition\n")
        for sample_id in sample_ids:
            condition = sample_id.split('_')[0]  # Extract condition from sample_id
            f.write(f"{sample_id},{condition}\n")
    
    # Create contrasts file
    with open(contrasts_file, 'w') as f:
        f.write("name,control,treatment\n")
        f.write("treatment_vs_control,control,treatment\n")

def create_readme(output_dir):
    """Create a README file for the test dataset"""
    readme_file = os.path.join(output_dir, 'README.md')
    
    with open(readme_file, 'w') as f:
        f.write("""# Synthetic RNA-Seq Pipeline Test Dataset

This directory contains a synthetic test dataset for the RNA-Seq Nextflow pipeline.

## Dataset Description

The test dataset consists of synthetic RNA-Seq data with paired-end reads.
The experiment compares gene expression between control and treatment conditions.

## Files

- **fastq/**: Contains paired-end FASTQ files for samples (control and treatment)
- **genome.fa**: Synthetic reference genome
- **annotation.gtf**: Synthetic gene annotation in GTF format
- **design.csv**: Sample metadata and experimental design
- **contrasts.csv**: Contrasts for differential expression analysis

## Usage

To run the pipeline with this test dataset:

```bash
nextflow run main.nf \\
  -profile test,docker \\
  --reads 'test_data_synthetic/fastq/*_R{1,2}.fastq.gz' \\
  --genome 'test_data_synthetic/genome.fa' \\
  --gtf 'test_data_synthetic/annotation.gtf' \\
  --design 'test_data_synthetic/design.csv' \\
  --contrasts 'test_data_synthetic/contrasts.csv' \\
  --outdir './results'
```

## Note

This is a synthetic dataset generated for testing purposes only.
For more realistic testing, use the real E. coli dataset in the `test_data` directory.
""")

def main():
    """Main function"""
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Generating synthetic test data in {args.output_dir}...")
    
    # Generate FASTQ files
    print("Generating FASTQ files...")
    sample_ids = generate_fastq_files(args.output_dir, args.num_samples, 
                                     args.num_reads, args.read_length)
    
    # Generate reference genome
    print("Generating reference genome...")
    genome_seq = generate_genome(args.output_dir, args.genome_size)
    
    # Generate gene annotation
    print("Generating gene annotation...")
    generate_annotation(args.output_dir, genome_seq, args.num_genes)
    
    # Generate design files
    print("Generating design files...")
    generate_design_files(args.output_dir, sample_ids)
    
    # Create README
    create_readme(args.output_dir)
    
    print("Synthetic test data generation complete!")

if __name__ == "__main__":
    main()
