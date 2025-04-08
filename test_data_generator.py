#!/usr/bin/env python3

import gzip
import random
import os

def generate_random_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def generate_quality_scores(length):
    return ''.join(chr(random.randint(33, 73)) for _ in range(length))

def create_fastq_file(filename, num_reads=1000, read_length=100):
    with gzip.open(filename, 'wt') as f:
        for i in range(num_reads):
            seq = generate_random_sequence(read_length)
            qual = generate_quality_scores(read_length)
            f.write(f'@READ_{i}\n')
            f.write(f'{seq}\n')
            f.write('+\n')
            f.write(f'{qual}\n')

def create_genome_file(filename, num_chromosomes=3, chr_length=10000):
    with open(filename, 'w') as f:
        for i in range(num_chromosomes):
            f.write(f'>chr{i+1}\n')
            sequence = generate_random_sequence(chr_length)
            # Write sequence in lines of 60 characters
            for j in range(0, len(sequence), 60):
                f.write(sequence[j:j+60] + '\n')

def create_gtf_file(filename, num_chromosomes=3, chr_length=10000, num_genes=50):
    with open(filename, 'w') as f:
        for i in range(num_genes):
            chr_num = random.randint(1, num_chromosomes)
            strand = random.choice(['+', '-'])
            gene_start = random.randint(1, chr_length - 1000)
            gene_end = gene_start + random.randint(500, 1000)
            
            # Write gene entry
            f.write(f'chr{chr_num}\ttest\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t'
                   f'gene_id "gene_{i}"; gene_name "gene_{i}";\n')
            
            # Add exons
            num_exons = random.randint(1, 3)
            current_pos = gene_start
            for j in range(num_exons):
                exon_start = current_pos
                exon_end = min(exon_start + random.randint(100, 300), gene_end)
                f.write(f'chr{chr_num}\ttest\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t'
                       f'gene_id "gene_{i}"; transcript_id "transcript_{i}_{j}"; exon_number "{j+1}";\n')
                current_pos = exon_end + random.randint(50, 100)
                if current_pos >= gene_end:
                    break

def main():
    os.makedirs('test_data', exist_ok=True)
    
    # Create paired-end FASTQ files
    create_fastq_file('test_data/sample1_R1.fastq.gz', num_reads=10000)
    create_fastq_file('test_data/sample1_R2.fastq.gz', num_reads=10000)
    
    # Create genome and annotation files
    create_genome_file('test_data/test_genome.fa')
    create_gtf_file('test_data/test_annotation.gtf')

if __name__ == '__main__':
    main() 