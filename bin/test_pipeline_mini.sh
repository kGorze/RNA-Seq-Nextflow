#!/bin/bash

# Exit on error
set -e

# Print commands
set -x

# Unset Python environment variables
unset PYTHONHOME
unset PYTHONPATH

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Change to project root
cd "$PROJECT_ROOT"

# Generate small test data
python3 test_data/test_data_generator.py --size small

# Run the complete pipeline with small test data
nextflow run main.nf \
    -profile test,docker \
    --reads 'test_data/sample*_R{1,2}.fastq.gz' \
    --genome 'test_data/test_genome.fa' \
    --gtf 'test_data/test_annotation.gtf' \
    --outdir './results' \
    --max_cpus 2 \
    --max_memory '4.GB'

# Return to original directory
cd - > /dev/null

echo "Small pipeline test completed successfully!"
