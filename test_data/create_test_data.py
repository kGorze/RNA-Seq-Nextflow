#!/usr/bin/env python3

import gzip
import os

def create_test_fastq(filename, num_reads=1000):
    """Create a test FASTQ file with realistic reads."""
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            # Create a realistic read header
            header = f"@INSTRUMENT:1:FLOWCELL:1:LANE:{i}:1:1:1:1"
            # Create a 100bp read with random ATCG
            seq = "".join(['ATCG'[i % 4] for i in range(100)])
            # Quality scores (Phred+33)
            qual = "".join(['F' for _ in range(100)])  # ASCII 70, good quality
            
            f.write(f"{header}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")

def create_test_genome(filename, size=1000):
    """Create a test genome FASTA file."""
    with open(filename, 'w') as f:
        f.write(">chr1\n")
        # Create a simple repeating sequence
        sequence = "".join(['ATCG'[i % 4] for i in range(size)])
        # Write sequence in lines of 60 characters
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

def create_test_gtf(filename):
    """Create a test GTF annotation file."""
    with open(filename, 'w') as f:
        f.write('##description: test GTF\n')
        # Add a few genes and transcripts
        for i in range(5):
            gene_start = i * 100 + 1
            gene_end = gene_start + 99
            f.write(f'chr1\ttest\tgene\t{gene_start}\t{gene_end}\t.\t+\t.\t'
                   f'gene_id "gene{i}"; gene_name "Gene{i}";\n')
            f.write(f'chr1\ttest\texon\t{gene_start}\t{gene_end}\t.\t+\t.\t'
                   f'gene_id "gene{i}"; transcript_id "trans{i}.1";\n')

def main():
    # Create test directory if it doesn't exist
    os.makedirs('test_data', exist_ok=True)
    
    # Create paired-end FASTQ files
    create_test_fastq('test_data/sample1_R1.fastq.gz')
    create_test_fastq('test_data/sample1_R2.fastq.gz')
    
    # Create genome and annotation files
    create_test_genome('test_data/test_genome.fa')
    create_test_gtf('test_data/test_annotation.gtf')

if __name__ == '__main__':
    main() 