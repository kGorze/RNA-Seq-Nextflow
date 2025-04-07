#!/bin/bash

# Script to create a final release package for the RNA-Seq Nextflow pipeline
# This script organizes all project files and creates a zip archive for delivery

# Set variables
PIPELINE_DIR=$(pwd)
RELEASE_DIR="${PIPELINE_DIR}/release"
RELEASE_FILE="rnaseq-nextflow-pipeline-v1.0.0.zip"

# Create release directory
mkdir -p "${RELEASE_DIR}"

# Create a manifest of all files
echo "Creating file manifest..."
find . -type f -not -path "*/\.*" -not -path "*/release/*" -not -path "*/results*/*" -not -path "*/work/*" -not -path "*/test_data*/*" | sort > "${RELEASE_DIR}/file_manifest.txt"

# Create a summary of the project
cat > "${RELEASE_DIR}/project_summary.md" << EOL
# RNA-Seq Nextflow Pipeline Project Summary

## Project Overview

This project implements a comprehensive RNA-Seq analysis pipeline using Nextflow. The pipeline provides end-to-end processing of RNA-Seq data from raw reads to differential expression results, with comprehensive reporting at each step.

## Key Features

- **Modular Design**: Each step is implemented as a separate module for flexibility and maintainability
- **Multiple Tool Options**: Support for different aligners (STAR, HISAT2) and quantification methods (Salmon, featureCounts)
- **Containerization**: Docker containers for reproducible execution across environments
- **Comprehensive Reporting**: Detailed HTML reports and visualizations of results
- **Scalability**: Efficient resource usage with Nextflow's execution engine
- **CI/CD Integration**: Automated testing and deployment with GitHub Actions

## Directory Structure

- **bin/**: Utility scripts for testing and data preparation
- **conf/**: Configuration files for different execution environments
- **docker/**: Dockerfiles and container build scripts
- **modules/**: Nextflow process modules for each pipeline component
- **.github/workflows/**: CI/CD workflow definitions
- **scripts/**: R scripts for analysis and visualization
- **test_data*/**: Test datasets for pipeline validation

## Getting Started

See the README.md file for detailed installation and usage instructions.

## License

This project is licensed under the MIT License.
EOL

# Create a version file
cat > "${RELEASE_DIR}/VERSION" << EOL
RNA-Seq Nextflow Pipeline
Version: 1.0.0
Release Date: $(date +"%Y-%m-%d")
EOL

# Create a zip archive of the project
echo "Creating release archive..."
zip -r "${RELEASE_DIR}/${RELEASE_FILE}" . \
    -x "*/\.*" \
    -x "*/release/*" \
    -x "*/results*/*" \
    -x "*/work/*" \
    -x "*/test_data*/*" \
    -x "*/\.*/*" \
    -x "*.git*" \
    -x "*.nextflow*"

# Create a minimal zip with just the essential files
echo "Creating minimal release archive..."
zip -r "${RELEASE_DIR}/rnaseq-nextflow-pipeline-minimal-v1.0.0.zip" \
    main.nf \
    nextflow.config \
    README.md \
    conf/ \
    modules/ \
    bin/test_pipeline_mini.sh \
    bin/generate_synthetic_data.py \
    docker/Dockerfile.*

echo "Release package created successfully!"
echo "Full release: ${RELEASE_DIR}/${RELEASE_FILE}"
echo "Minimal release: ${RELEASE_DIR}/rnaseq-nextflow-pipeline-minimal-v1.0.0.zip"
echo "Project summary: ${RELEASE_DIR}/project_summary.md"
