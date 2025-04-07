#!/bin/bash

# Script to test the RNA-Seq Nextflow pipeline with the mini test dataset
# This script runs the pipeline with the minimal test dataset for quick validation

# Set variables
PIPELINE_DIR=$(pwd)
TEST_DATA_DIR="${PIPELINE_DIR}/test_data_mini"
RESULTS_DIR="${PIPELINE_DIR}/results_mini_test"

# Check if test data exists, if not prepare it
if [ ! -d "${TEST_DATA_DIR}" ]; then
    echo "Test data directory not found. Preparing mini test dataset..."
    bash bin/prepare_mini_test_data.sh
fi

# Run the pipeline with the mini test dataset
echo "Running RNA-Seq pipeline with mini test dataset..."
nextflow run main.nf \
  -profile test,docker \
  --reads "${TEST_DATA_DIR}/fastq/*_R{1,2}.fastq.gz" \
  --genome "${TEST_DATA_DIR}/genome.fa" \
  --gtf "${TEST_DATA_DIR}/annotation.gtf" \
  --design "${TEST_DATA_DIR}/design.csv" \
  --contrasts "${TEST_DATA_DIR}/contrasts.csv" \
  --outdir "${RESULTS_DIR}" \
  --max_cpus 2 \
  --max_memory '4.GB'

# Check if the pipeline completed successfully
if [ $? -eq 0 ]; then
    echo "Pipeline test with mini dataset completed successfully!"
    
    # Check for key output files
    echo "Checking for key output files..."
    
    # Define expected output directories
    expected_dirs=(
        "fastqc"
        "multiqc"
        "trimmed"
        "alignment"
        "quantification"
        "deseq2"
        "final_report"
    )
    
    # Check if each expected directory exists
    all_found=true
    for dir in "${expected_dirs[@]}"; do
        if [ ! -d "${RESULTS_DIR}/${dir}" ]; then
            echo "ERROR: Expected output directory ${dir} not found!"
            all_found=false
        else
            echo "Found output directory: ${dir}"
        fi
    done
    
    if $all_found; then
        echo "All expected output directories found. Test passed!"
    else
        echo "Some expected output directories are missing. Test failed!"
        exit 1
    fi
else
    echo "Pipeline test with mini dataset failed!"
    exit 1
fi
