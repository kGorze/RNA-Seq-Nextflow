#!/usr/bin/env nextflow

/*
 * Test data generation module for RNA-Seq pipeline
 * Creates synthetic test data for pipeline validation
 */

process GENERATE_TEST_DATA {
    label 'process_low'
    
    publishDir "${params.outdir}/test_data", mode: 'copy'
    
    output:
    path "sample*_R{1,2}.fastq.gz", emit: reads
    path "test_genome.fa", emit: genome
    path "test_annotation.gtf", emit: gtf
    
    script:
    """
    #!/usr/bin/env python3
    import gzip
    import random
    
    # Function to generate random DNA sequence
    def random_dna(length):
        return ''.join(random.choice('ACGT') for _ in range(length))
    
    # Function to generate FASTQ entry
    def generate_fastq_entry(read_id, seq, qual_score=30):
        quality = chr(33 + qual_score) * len(seq)
        return f"@{read_id}\\n{seq}\\n+\\n{quality}\\n"
    
    # Generate paired-end reads for two samples
    for sample_id in [1, 2]:
        # Create R1 file
        with gzip.open(f"sample{sample_id}_R1.fastq.gz", 'wt') as r1:
            for i in range(1000):  # 1000 reads per sample
                seq = random_dna(100)
                r1.write(generate_fastq_entry(f"READ_{sample_id}_{i}", seq))
        
        # Create R2 file
        with gzip.open(f"sample{sample_id}_R2.fastq.gz", 'wt') as r2:
            for i in range(1000):  # 1000 reads per sample
                seq = random_dna(100)
                r2.write(generate_fastq_entry(f"READ_{sample_id}_{i}", seq))
    
    # Generate a small test genome
    with open("test_genome.fa", 'w') as genome:
        # Create 3 chromosomes
        for chrom in range(1, 4):
            genome.write(f">chr{chrom}\\n")
            # Each chromosome is 10kb
            seq = random_dna(10000)
            # Write sequence in lines of 80 characters
            for i in range(0, len(seq), 80):
                genome.write(seq[i:i+80] + "\\n")
    
    # Generate a test GTF annotation
    with open("test_annotation.gtf", 'w') as gtf:
        # Header
        gtf.write('##description: Test annotation for RNA-Seq pipeline\\n')
        gtf.write('#!genome-build Test\\n')
        
        # Create 10 genes with 2 transcripts each
        for gene_id in range(1, 11):
            chrom = random.randint(1, 3)
            strand = random.choice(['+', '-'])
            gene_start = random.randint(1000, 9000)
            gene_end = gene_start + random.randint(500, 1000)
            
            # Gene entry
            gtf.write(f'chr{chrom}\\ttest\\tgene\\t{gene_start}\\t{gene_end}\\t.\\t{strand}\\t.\\t')
            gtf.write(f'gene_id "gene{gene_id}"; gene_name "Gene{gene_id}";\\n')
            
            # Create 2 transcripts for each gene
            for transcript_id in range(1, 3):
                tx_id = f"{gene_id}_{transcript_id}"
                tx_start = gene_start
                tx_end = gene_end
                
                # Transcript entry
                gtf.write(f'chr{chrom}\\ttest\\ttranscript\\t{tx_start}\\t{tx_end}\\t.\\t{strand}\\t.\\t')
                gtf.write(f'gene_id "gene{gene_id}"; transcript_id "tx{tx_id}"; gene_name "Gene{gene_id}";\\n')
                
                # Create 3-5 exons for each transcript
                num_exons = random.randint(3, 5)
                exon_start = tx_start
                
                for exon_id in range(1, num_exons + 1):
                    exon_length = random.randint(100, 200)
                    exon_end = exon_start + exon_length
                    
                    if exon_end > tx_end:
                        exon_end = tx_end
                    
                    # Exon entry
                    gtf.write(f'chr{chrom}\\ttest\\texon\\t{exon_start}\\t{exon_end}\\t.\\t{strand}\\t.\\t')
                    gtf.write(f'gene_id "gene{gene_id}"; transcript_id "tx{tx_id}"; exon_id "exon{tx_id}_{exon_id}"; gene_name "Gene{gene_id}";\\n')
                    
                    exon_start = exon_end + random.randint(50, 100)
                    if exon_start >= tx_end:
                        break
    """
}
