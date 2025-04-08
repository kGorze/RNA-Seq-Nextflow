#!/usr/bin/env python3

import gzip
import random
import os
import numpy as np
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_random_sequence(length):
    """Generate a random DNA sequence with realistic GC content."""
    gc_content = random.uniform(0.4, 0.6)  # Realistic GC content range
    gc_bases = int(length * gc_content)
    at_bases = length - gc_bases
    
    sequence = ['G', 'C'] * (gc_bases // 2) + ['A', 'T'] * (at_bases // 2)
    if gc_bases % 2:
        sequence.append(random.choice(['G', 'C']))
    if at_bases % 2:
        sequence.append(random.choice(['A', 'T']))
    
    random.shuffle(sequence)
    return ''.join(sequence)

def generate_quality_scores(length, mean_quality=35, std_quality=5):
    """Generate realistic quality scores following a normal distribution."""
    qualities = np.random.normal(mean_quality, std_quality, length)
    qualities = np.clip(qualities, 0, 40)  # Limit quality scores between 0 and 40
    return ''.join(chr(int(q) + 33) for q in qualities)  # Convert to ASCII

def create_fastq_file(filename, num_reads=10000, read_length=100, paired_end=True):
    """Create realistic FASTQ files with proper read structure and quality scores."""
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            # Generate sequence with realistic GC content
            seq = generate_random_sequence(read_length)
            qual = generate_quality_scores(read_length)
            
            # Write FASTQ entry
            f.write(f'@READ_{i}\n')
            f.write(f'{seq}\n')
            f.write('+\n')
            f.write(f'{qual}\n')

def create_genome_file(filename, num_chromosomes=3, chr_length=100000):
    """Create a genome file with realistic chromosome structure."""
    with open(filename, 'w') as f:
        for i in range(num_chromosomes):
            f.write(f'>chr{i+1}\n')
            sequence = generate_random_sequence(chr_length)
            # Write sequence in lines of 60 characters
            for j in range(0, len(sequence), 60):
                f.write(sequence[j:j+60] + '\n')

def create_gtf_file(filename, num_chromosomes=3, chr_length=100000, num_genes=50):
    """Create a GTF file with realistic gene structure including UTRs and exons."""
    with open(filename, 'w') as f:
        for i in range(num_genes):
            chr_num = random.randint(1, num_chromosomes)
            strand = random.choice(['+', '-'])
            
            # Generate gene coordinates ensuring they fit within chromosome
            gene_start = random.randint(1000, chr_length - 5000)
            gene_length = random.randint(1000, 5000)
            gene_end = min(gene_start + gene_length, chr_length - 1000)
            
            # Write gene entry
            f.write(f'chr{chr_num}\ttest\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t'
                   f'gene_id "gene_{i}"; gene_name "gene_{i}"; gene_biotype "protein_coding";\n')
            
            # Add 5' UTR
            utr5_start = gene_start
            utr5_end = gene_start + random.randint(50, 200)
            f.write(f'chr{chr_num}\ttest\tfive_prime_utr\t{utr5_start}\t{utr5_end}\t.\t{strand}\t.\t'
                   f'gene_id "gene_{i}"; transcript_id "transcript_{i}";\n')
            
            # Add exons
            num_exons = random.randint(2, 5)
            current_pos = utr5_end + 1
            exon_number = 1
            
            for j in range(num_exons):
                exon_start = current_pos
                exon_length = random.randint(100, 300)
                exon_end = min(exon_start + exon_length, gene_end - 200)  # Leave space for 3' UTR
                
                f.write(f'chr{chr_num}\ttest\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t'
                       f'gene_id "gene_{i}"; transcript_id "transcript_{i}"; exon_number "{exon_number}";\n')
                
                if j < num_exons - 1:  # Add intron if not the last exon
                    intron_start = exon_end + 1
                    intron_length = random.randint(50, 500)
                    current_pos = intron_start + intron_length
                else:
                    current_pos = exon_end + 1
                
                exon_number += 1
            
            # Add 3' UTR
            utr3_start = current_pos
            utr3_end = gene_end
            f.write(f'chr{chr_num}\ttest\tthree_prime_utr\t{utr3_start}\t{utr3_end}\t.\t{strand}\t.\t'
                   f'gene_id "gene_{i}"; transcript_id "transcript_{i}";\n')

def main():
    parser = argparse.ArgumentParser(description='Generate test data files for bioinformatics analysis')
    parser.add_argument('--only-fastq', action='store_true', help='Generate only FASTQ files')
    parser.add_argument('--only-genome', action='store_true', help='Generate only genome file')
    parser.add_argument('--only-annotation', action='store_true', help='Generate only annotation file')
    args = parser.parse_args()

    # If no specific option is selected, generate all files
    if not any([args.only_fastq, args.only_genome, args.only_annotation]):
        args.only_fastq = True
        args.only_genome = True
        args.only_annotation = True

    generated_files = []

    if args.only_fastq:
        # Create paired-end FASTQ files with realistic parameters
        create_fastq_file('sample1_R1.fastq.gz', num_reads=10000, read_length=100)
        create_fastq_file('sample1_R2.fastq.gz', num_reads=10000, read_length=100)
        generated_files.extend(['sample1_R1.fastq.gz', 'sample1_R2.fastq.gz'])
    
    if args.only_genome:
        # Create genome file with realistic structure
        create_genome_file('test_genome.fa', num_chromosomes=3, chr_length=100000)
        generated_files.append('test_genome.fa')
    
    if args.only_annotation:
        # Create annotation file with realistic structure
        create_gtf_file('test_annotation.gtf', num_chromosomes=3, chr_length=100000, num_genes=50)
        generated_files.append('test_annotation.gtf')
    
    print("Test data generation completed successfully!")
    print("Generated files:")
    for file in generated_files:
        print(f"- {file}")

if __name__ == '__main__':
    main() 