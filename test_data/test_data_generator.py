#!/usr/bin/env python3

import gzip
import random
import os
import numpy as np
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Size configurations
SIZE_CONFIGS = {
    "tiny": {
        "num_reads": 200,
        "read_length": 50,
        "num_chromosomes": 1,
        "chr_length": 2000,
        "num_genes": 2,
    },
    "mini": {
        "num_reads": 1000,
        "read_length": 50,
        "num_chromosomes": 1,
        "chr_length": 10000,
        "num_genes": 10,
    },
    "small": {
        "num_reads": 2000,
        "read_length": 50,
        "num_chromosomes": 1,
        "chr_length": 20000,
        "num_genes": 15,
    },
    "medium": {
        "num_reads": 10000,
        "read_length": 100,
        "num_chromosomes": 3,
        "chr_length": 100000,
        "num_genes": 50,
    },
    "large": {
        "num_reads": 20000,
        "read_length": 150,
        "num_chromosomes": 5,
        "chr_length": 200000,
        "num_genes": 100,
    },
}


def generate_random_sequence(length):
    """Generate a random DNA sequence with realistic GC content."""
    gc_content = random.uniform(0.4, 0.6)  # Realistic GC content range
    gc_bases = int(length * gc_content)
    at_bases = length - gc_bases

    sequence = ["G", "C"] * (gc_bases // 2) + ["A", "T"] * (at_bases // 2)
    if gc_bases % 2:
        sequence.append(random.choice(["G", "C"]))
    if at_bases % 2:
        sequence.append(random.choice(["A", "T"]))

    random.shuffle(sequence)
    return "".join(sequence)


def generate_quality_scores(length, mean_quality=35, std_quality=5):
    """Generate realistic quality scores following a normal distribution."""
    qualities = np.random.normal(mean_quality, std_quality, length)
    qualities = np.clip(qualities, 0, 40)  # Limit quality scores between 0 and 40
    return "".join(chr(int(q) + 33) for q in qualities)  # Convert to ASCII


def create_fastq_file(filename, size="medium"):
    """Create realistic FASTQ files with proper read structure and quality scores."""
    config = SIZE_CONFIGS[size]
    with gzip.open(filename, "wt") as f:
        for i in range(config["num_reads"]):
            seq = generate_random_sequence(config["read_length"])
            qual = generate_quality_scores(config["read_length"])
            f.write(f"@READ_{i}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")


def create_genome_file(filename, size="medium"):
    """Create a genome file with realistic chromosome structure."""
    config = SIZE_CONFIGS[size]
    with open(filename, "w") as f:
        for i in range(config["num_chromosomes"]):
            f.write(f">chr{i+1}\n")
            sequence = generate_random_sequence(config["chr_length"])
            for j in range(0, len(sequence), 60):
                f.write(sequence[j : j + 60] + "\n")


def create_gtf_file(filename, size="medium"):
    """Create a GTF file with realistic gene structure including transcripts, UTRs and exons."""
    config = SIZE_CONFIGS[size]
    with open(filename, "w") as f:
        for i in range(config["num_genes"]):
            chr_num = random.randint(1, config["num_chromosomes"])
            strand = random.choice(["+", "-"])

            # Calculate gene coordinates with more space
            if size == "tiny":
                # For tiny dataset, use smaller ranges
                gene_start = random.randint(
                    100, config["chr_length"] - 500
                )  # Start further in
                gene_length = random.randint(200, 400)  # Smaller minimum gene length
                gene_end = min(gene_start + gene_length, config["chr_length"] - 100)
            else:
                # Original logic for larger datasets
                gene_start = random.randint(1000, config["chr_length"] - 15000)
                gene_length = random.randint(5000, 10000)
                gene_end = min(gene_start + gene_length, config["chr_length"] - 1000)

            # Write gene entry
            f.write(
                f"chr{chr_num}\ttest\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
                f'gene_id "gene_{i}"; gene_name "gene_{i}"; gene_biotype "protein_coding";\n'
            )

            # Write transcript entry (before its subfeatures)
            f.write(
                f"chr{chr_num}\ttest\ttranscript\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
                f'gene_id "gene_{i}"; transcript_id "transcript_{i}"; gene_name "gene_{i}";\n'
            )

            # Calculate feature lengths
            total_length = gene_end - gene_start + 1
            if size == "tiny":
                utr_length = total_length // 5  # 20% for each UTR in tiny dataset
            else:
                utr_length = total_length // 10  # 10% for each UTR in larger datasets
            exon_region_length = total_length - (
                2 * utr_length
            )  # Remaining space for exons

            # Add 5' UTR
            utr5_start = gene_start
            utr5_end = utr5_start + utr_length - 1
            f.write(
                f"chr{chr_num}\ttest\tfive_prime_utr\t{utr5_start}\t{utr5_end}\t.\t{strand}\t.\t"
                f'gene_id "gene_{i}"; transcript_id "transcript_{i}";\n'
            )

            # Add exons
            if size == "tiny":
                num_exons = random.randint(2, 3)  # Fewer exons for tiny dataset
            else:
                num_exons = random.randint(3, 5)
            exon_length = exon_region_length // num_exons
            intron_length = (exon_region_length - (exon_length * num_exons)) // (
                num_exons - 1
            )

            current_pos = utr5_end + 1
            for j in range(num_exons):
                exon_start = current_pos
                exon_end = exon_start + exon_length - 1

                f.write(
                    f"chr{chr_num}\ttest\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t"
                    f'gene_id "gene_{i}"; transcript_id "transcript_{i}"; exon_number "{j+1}";\n'
                )

                if j < num_exons - 1:
                    current_pos = exon_end + intron_length + 1
                else:
                    current_pos = exon_end + 1

            # Add 3' UTR
            utr3_start = current_pos
            utr3_end = gene_end
            f.write(
                f"chr{chr_num}\ttest\tthree_prime_utr\t{utr3_start}\t{utr3_end}\t.\t{strand}\t.\t"
                f'gene_id "gene_{i}"; transcript_id "transcript_{i}";\n'
            )


def main():
    parser = argparse.ArgumentParser(
        description="Generate test data for RNA-Seq pipeline"
    )
    parser.add_argument(
        "--only-fastq", action="store_true", help="Generate only FASTQ files"
    )
    parser.add_argument(
        "--only-genome", action="store_true", help="Generate only genome file"
    )
    parser.add_argument(
        "--only-annotation", action="store_true", help="Generate only annotation file"
    )
    parser.add_argument(
        "--size",
        choices=["tiny", "mini", "small", "medium", "large"],
        default="medium",
        help="Size of the test dataset (default: medium)",
    )
    args = parser.parse_args()

    # If no specific option is selected, generate all files
    if not any([args.only_fastq, args.only_genome, args.only_annotation]):
        args.only_fastq = True
        args.only_genome = True
        args.only_annotation = True

    print(f"Generating {args.size} test dataset...")

    if args.only_fastq:
        create_fastq_file(f"test_data/sample1_R1.fastq.gz", size=args.size)
        create_fastq_file(f"test_data/sample1_R2.fastq.gz", size=args.size)
        create_fastq_file(f"test_data/sample2_R1.fastq.gz", size=args.size)
        create_fastq_file(f"test_data/sample2_R2.fastq.gz", size=args.size)

    if args.only_genome:
        create_genome_file(f"test_data/test_genome.fa", size=args.size)

    if args.only_annotation:
        create_gtf_file(f"test_data/test_annotation.gtf", size=args.size)

    # Create design and contrasts files
    if args.only_fastq:  # Only create these if we're generating FASTQ files
        with open("test_data/design.csv", "w") as f:
            f.write("sample_id,condition\n")
            f.write("sample1,control\n")
            f.write("sample2,treatment\n")

        with open("test_data/contrasts.csv", "w") as f:
            f.write("name,control,treatment\n")
            f.write("treatment_vs_control,control,treatment\n")

    print(f"Test dataset generation complete!")
    config = SIZE_CONFIGS[args.size]
    print(f"\nDataset details:")
    print(f"- Number of reads per file: {config['num_reads']}")
    print(f"- Read length: {config['read_length']} bp")
    print(f"- Number of chromosomes: {config['num_chromosomes']}")
    print(f"- Chromosome length: {config['chr_length']} bp")
    print(f"- Number of genes: {config['num_genes']}")


if __name__ == "__main__":
    main()
