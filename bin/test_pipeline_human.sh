#!/bin/bash

# Script to test the RNA-Seq Nextflow pipeline with the human test dataset
# This script runs the pipeline with the human test dataset (chromosome 22) for validation

# Set variables
PIPELINE_DIR=$(pwd)
TEST_DATA_DIR="${PIPELINE_DIR}/test_data_human"
RESULTS_DIR="${PIPELINE_DIR}/results_human_test"

# Check if test data exists, if not prepare it
if [ ! -d "${TEST_DATA_DIR}" ]; then
    echo "Test data directory not found. Preparing human test dataset..."
    bash bin/prepare_human_test_data.sh
fi

# Run the pipeline with the human test dataset
echo "Running RNA-Seq pipeline with human test dataset..."
nextflow run main.nf \
  -profile test,docker \
  --reads "${TEST_DATA_DIR}/fastq/*_R{1,2}.fastq.gz" \
  --genome "${TEST_DATA_DIR}/genome.fa" \
  --gtf "${TEST_DATA_DIR}/annotation.gtf" \
  --design "${TEST_DATA_DIR}/design.csv" \
  --contrasts "${TEST_DATA_DIR}/contrasts.csv" \
  --outdir "${RESULTS_DIR}" \
  --max_cpus 8 \
  --max_memory '16.GB'

# Check if the pipeline completed successfully
if [ $? -eq 0 ]; then
    echo "Pipeline test with human dataset completed successfully!"
    
    # Check for key output files
    echo "Checking for key output files..."
    
    # Define expected output directories and files
    expected_dirs=(
        "fastqc"
        "multiqc"
        "trimmed"
        "alignment"
        "quantification"
        "deseq2"
        "final_report"
    )
    
    expected_files=(
        "final_report/final_report.html"
        "deseq2/results/treatment_vs_control_results.csv"
        "multiqc/multiqc_report.html"
        "quantification/featurecounts/merged/merged_gene_counts.csv"
        "deseq2/report/deseq2_report.html"
    )
    
    # Check if each expected directory exists
    all_dirs_found=true
    for dir in "${expected_dirs[@]}"; do
        if [ ! -d "${RESULTS_DIR}/${dir}" ]; then
            echo "ERROR: Expected output directory ${dir} not found!"
            all_dirs_found=false
        else
            echo "Found output directory: ${dir}"
        fi
    done
    
    # Check if each expected file exists
    all_files_found=true
    for file in "${expected_files[@]}"; do
        if [ ! -f "${RESULTS_DIR}/${file}" ]; then
            echo "ERROR: Expected output file ${file} not found!"
            all_files_found=false
        else
            echo "Found output file: ${file}"
        fi
    done
    
    if $all_dirs_found && $all_files_found; then
        echo "All expected output directories and files found. Test passed!"
    else
        echo "Some expected outputs are missing. Test failed!"
        exit 1
    fi
else
    echo "Pipeline test with human dataset failed!"
    exit 1
fi
