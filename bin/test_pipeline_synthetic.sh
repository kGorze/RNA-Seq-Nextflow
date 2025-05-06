#!/bin/bash

set -e  # Exit on error
set -u  # Exit on undefined variable

echo "Running synthetic data pipeline test..."

# Unset Python environment variables
unset PYTHONHOME
unset PYTHONPATH

# Generate synthetic test data
echo "Generating synthetic test dataset..."
python test_data/test_data_generator.py --only-fastq
python test_data/test_data_generator.py --only-genome
python test_data/test_data_generator.py --only-annotation

# Run pipeline with synthetic dataset
echo "Running pipeline with synthetic dataset..."
nextflow run main.nf \
    -profile test,docker \
    --reads 'test_data/sample*_R{1,2}.fastq.gz' \
    --genome 'test_data/test_genome.fa' \
    --gtf 'test_data/test_annotation.gtf' \
    --design 'test_data/design.csv' \
    --contrasts 'test_data/contrasts.csv' \
    --outdir './results' \
    --max_cpus 2 \
    --max_memory '4.GB'

echo "Synthetic data test completed successfully!"
