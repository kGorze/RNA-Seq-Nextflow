#!/bin/bash

# Build script for RNA-Seq pipeline Docker containers
# This script builds all the Docker containers required for the pipeline

# Set variables
DOCKER_REPO="rnaseq-pipeline"
VERSION="1.0.0"

# Create adapter files for Trimmomatic
echo "Creating adapter files for Trimmomatic..."
cat > TruSeq3-PE.fa << EOL
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
EOL

cat > TruSeq3-SE.fa << EOL
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
EOL

# Build base image
echo "Building base image..."
docker build -t ${DOCKER_REPO}/base:${VERSION} -f Dockerfile.base .

# Build tool-specific images
echo "Building FastQC image..."
docker build -t ${DOCKER_REPO}/fastqc:${VERSION} -f Dockerfile.fastqc .

echo "Building MultiQC image..."
docker build -t ${DOCKER_REPO}/multiqc:${VERSION} -f Dockerfile.multiqc .

echo "Building Trimmomatic image..."
docker build -t ${DOCKER_REPO}/trimmomatic:${VERSION} -f Dockerfile.trimmomatic .

echo "Building STAR image..."
docker build -t ${DOCKER_REPO}/star:${VERSION} -f Dockerfile.star .

echo "Building HISAT2 image..."
docker build -t ${DOCKER_REPO}/hisat2:${VERSION} -f Dockerfile.hisat2 .

echo "Building Salmon image..."
docker build -t ${DOCKER_REPO}/salmon:${VERSION} -f Dockerfile.salmon .

echo "Building featureCounts image..."
docker build -t ${DOCKER_REPO}/featurecounts:${VERSION} -f Dockerfile.featurecounts .

echo "Building DESeq2 image..."
docker build -t ${DOCKER_REPO}/deseq2:${VERSION} -f Dockerfile.deseq2 .

echo "All Docker images built successfully!"
echo "To push images to Docker Hub, run:"
echo "docker login"
echo "docker push ${DOCKER_REPO}/base:${VERSION}"
echo "docker push ${DOCKER_REPO}/fastqc:${VERSION}"
echo "docker push ${DOCKER_REPO}/multiqc:${VERSION}"
echo "docker push ${DOCKER_REPO}/trimmomatic:${VERSION}"
echo "docker push ${DOCKER_REPO}/star:${VERSION}"
echo "docker push ${DOCKER_REPO}/hisat2:${VERSION}"
echo "docker push ${DOCKER_REPO}/salmon:${VERSION}"
echo "docker push ${DOCKER_REPO}/featurecounts:${VERSION}"
echo "docker push ${DOCKER_REPO}/deseq2:${VERSION}"
