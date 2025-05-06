# RNA-Seq Nextflow Pipeline v1.0.0

[![CI](https://github.com/username/rnaseq-nextflow-pipeline/workflows/CI/badge.svg)](https://github.com/username/rnaseq-nextflow-pipeline/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/docker/automated/username/rnaseq-pipeline.svg)](https://hub.docker.com/r/username/rnaseq-pipeline)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A reproducible, scalable, and containerized RNA-Seq analysis pipeline implemented in Nextflow.

## Overview

This pipeline performs RNA-Seq analysis from raw reads to differential expression results. It includes quality control, read trimming, alignment, quantification, and differential expression analysis, with comprehensive reporting at each step.

### Pipeline Features

- **Modular Design**: Each step is implemented as a separate module for flexibility and maintainability
- **Multiple Tool Options**: Support for different aligners (STAR, HISAT2) and quantification methods (Salmon, featureCounts)
- **Containerization**: Docker containers for reproducible execution across environments
- **Comprehensive Reporting**: Detailed HTML reports and visualizations of results
- **Scalability**: Efficient resource usage with Nextflow's execution engine
- **CI/CD Integration**: Automated testing and deployment with GitHub Actions

## Quick Start

```bash
# Clone the repository
git clone https://github.com/kGorze/RNA-Seq-Nextflow
cd RNA-Seq-Nextflow

# Run the pipeline with test data
nextflow run main.nf -profile test,docker

# Run with your own data
nextflow run main.nf \
  -profile docker \
  --reads '/path/to/reads/*_R{1,2}.fastq.gz' \
  --genome '/path/to/genome.fa' \
  --gtf '/path/to/annotation.gtf' \
  --design '/path/to/design.csv' \
  --contrasts '/path/to/contrasts.csv' \
  --outdir './results'
```

## Installation

### Requirements

- [Nextflow](https://www.nextflow.io/) (>=20.10.0)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) (for containerized execution)
- Java 8 or later

### Installation Steps

1. Install Nextflow:
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. Add Nextflow to your PATH:
   ```bash
   mv nextflow /usr/local/bin/
   ```

3. Install Docker:
   ```bash
   curl -fsSL https://get.docker.com | sh
   ```

4. Clone the repository:
   ```bash
   git clone https://github.com/kGorze/RNA-Seq-Nextflow.git
   cd rnaseq-nextflow-pipeline
   ```

## Pipeline Usage

### Basic Usage

```bash
nextflow run main.nf \
  -profile docker \
  --reads '/path/to/reads/*_R{1,2}.fastq.gz' \
  --genome '/path/to/genome.fa' \
  --gtf '/path/to/annotation.gtf' \
  --outdir './results'
```

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--reads` | Path to input FASTQ files (must be quoted and include the {1,2} pattern for paired-end data) |
| `--genome` | Path to reference genome FASTA file |
| `--gtf` | Path to gene annotation GTF file |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory path | `./results` |
| `--design` | Path to experimental design CSV file | `null` |
| `--contrasts` | Path to contrasts CSV file for differential expression | `null` |
| `--aligner` | Alignment tool to use (star, hisat2) | `star` |
| `--quantification` | Quantification method to use (salmon, featurecounts) | `featurecounts` |
| `--skip_trimming` | Skip the trimming step | `false` |
| `--skip_alignment` | Skip the alignment step (only valid with salmon) | `false` |
| `--skip_qc` | Skip the quality control step | `false` |
| `--skip_deseq2` | Skip the differential expression analysis | `false` |
| `--max_cpus` | Maximum number of CPUs to use | `16` |
| `--max_memory` | Maximum memory to use | `128.GB` |

### Profiles

The pipeline comes with several configuration profiles:

- `docker`: Uses Docker containers for execution
- `singularity`: Uses Singularity containers for execution
- `test`: Runs the pipeline with the test dataset
- `standard`: Runs the pipeline without containers (requires all tools to be installed)
- `slurm`: Adds configuration for Slurm workload manager
- `sge`: Adds configuration for Sun Grid Engine

Profiles can be combined using comma-separated values:
```bash
nextflow run main.nf -profile test,docker
```

### Commands guide-book

1. First we've got the pipeline execution commands:
```
# Run the pipeline with default parameters
nextflow run main.nf

# Run with test profile and Docker
nextflow run main.nf -profile test,docker

# Resume a previous run (recommended)
nextflow run main.nf -profile test,docker -resume

# Run with custom configuration
nextflow run main.nf -c custom.config

# Run with specific work directory
nextflow run main.nf -w /path/to/work/dir

# Run pipeline with debug log level
nextflow run main.nf -profile test,docker -log debug

# View pipeline help
nextflow run main.nf --help

# View pipeline version and parameters
nextflow run main.nf --version
```


2. Then, when we want to change the default configuration we can use those commands:
```
# Run with different profiles
nextflow run main.nf -profile docker
nextflow run main.nf -profile singularity
nextflow run main.nf -profile conda

# Specify custom resources
nextflow run main.nf --max_memory '64.GB' --max_cpus 16

# Override parameters
nextflow run main.nf --reads '/path/to/reads/*_{1,2}.fastq.gz'
nextflow run main.nf --genome '/path/to/genome.fa'
nextflow run main.nf --annotation '/path/to/annotation.gtf'
```

3. If you don't have the specified data from your experiment you can use for learning purposes the synthetic one which we provide with those commands:

```
# Generate test data
cd test_data
python test_data_generator.py

# Clean test data
rm test_data/*.fastq.gz test_data/*.fa test_data/*.gtf

# Regenerate specific test files
python test_data_generator.py --only-fastq  # (#TODO)
python test_data_generator.py --only-genome  # (#TODO)
python test_data_generator.py --only-annotation  # (#TODO)
```

4. After those initial settings we can continue with the managements commands:

```
# Clean work directory
rm -rf work/

# Clean results directory
rm -rf results/

# View pipeline execution report(here you can use the browser of your choice, all of those html are based on the dependencies which are not made by me)
firefox results/multiqc/multiqc_report.html

# List all processes in the work directory
ls -l work/*/

# View Nextflow logs
cat .nextflow.log

# View specific process logs
cat work/XX/XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/.command.log
```

5. Docker/Container management commands:
```
# Pull required Docker images
docker pull rocker/r-base:4.1.0
docker pull biocontainers/star:2.7.10a--h9ee0642_0
docker pull biocontainers/fastqc:0.11.9--0
docker pull biocontainers/trimmomatic:0.39--hdfd78af_2
docker pull biocontainers/subread:2.0.1--hed695b0_0

# List downloaded images
docker images

# Remove unused images
docker system prune
```

6. We also have some debug commands, if anything would brake:
```
# Run pipeline with debug log level
nextflow run main.nf -profile test,docker -log debug

# Inspect a specific process
cd work/XX/XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/
bash .command.run

# View process environment
cat .command.env

# View process script
cat .command.sh
```
7. Maintain all the information about the performance:
```
# Monitor resource usage
top -u $USER

# Monitor disk usage
du -sh work/ results/

# View Nextflow runtime statistics
nextflow log
```
8. Eventually we can watch the results with our file viewer:
```
# View merged gene counts
cat results/featurecounts/merged/merged_gene_counts.csv

# View merged transcript counts
cat results/featurecounts/merged/merged_transcript_counts.csv

# View alignment statistics
cat results/star/*/Log.final.out

# View FastQC reports
firefox results/fastqc/*.html
```

## Input Files

### Read Files

The pipeline accepts paired-end FASTQ files. Files should be specified using the `--reads` parameter with a glob pattern:

```bash
--reads '/path/to/reads/*_R{1,2}.fastq.gz'
```

This pattern assumes that paired files have the same prefix and are identified by `_R1` and `_R2` suffixes.

### Experimental Design

For differential expression analysis, you need to provide a design file in CSV format with at least the following columns:

```csv
sample_id,condition
sample1,control
sample2,control
sample3,treatment
sample4,treatment
```

- `sample_id`: Must match the prefix of your FASTQ files
- `condition`: Experimental condition for each sample

Additional columns can be included for more complex designs (e.g., batch, sex, age).

### Contrasts

To specify which conditions to compare in the differential expression analysis, provide a contrasts file in CSV format:

```csv
name,control,treatment
treatment_vs_control,control,treatment
```

- `name`: Name for the contrast (used in output files)
- `control`: Reference condition
- `treatment`: Test condition

## Output Files

The pipeline organizes outputs into the following directory structure:

```
results/
├── 01_quality_control/
│   ├── fastqc/
│   └── multiqc/
├── 02_trimming/
├── 03_alignment/
│   ├── star/ or hisat2/
│   └── stats/
├── 04_quantification/
│   ├── featurecounts/ or salmon/
│   └── merged/
├── 05_differential_expression/
│   ├── deseq2/
│   └── visualizations/
├── 06_reports/
│   └── final_report.html
└── pipeline_info/
```

### Key Output Files

| File | Description |
|------|-------------|
| `01_quality_control/multiqc/multiqc_report.html` | Quality control summary report |
| `03_alignment/stats/alignment_stats.csv` | Summary of alignment statistics |
| `04_quantification/merged/merged_gene_counts.csv` | Merged gene count matrix |
| `05_differential_expression/deseq2/results/*_results.csv` | Differential expression results |
| `06_reports/final_report.html` | Comprehensive analysis report |

## Test Datasets

The pipeline includes several test datasets for validation and testing:

1. **Minimal Test Dataset**: Very small dataset for quick CI/CD testing
   ```bash
   bash bin/prepare_mini_test_data.sh
   nextflow run main.nf -profile test,docker --reads 'test_data_mini/fastq/*_R{1,2}.fastq.gz'
   ```

2. **E. coli Test Dataset**: Small bacterial dataset for basic validation
   ```bash
   bash bin/prepare_test_data.sh
   nextflow run main.nf -profile test,docker --reads 'test_data/fastq/*_R{1,2}.fastq.gz'
   ```

3. **Synthetic Test Dataset**: Generated data with controllable parameters
   ```bash
   python bin/generate_synthetic_data.py
   nextflow run main.nf -profile test,docker --reads 'test_data_synthetic/fastq/*_R{1,2}.fastq.gz'
   ```

4. **Human Test Dataset**: Subset of human RNA-Seq data (chromosome 22)
   ```bash
   bash bin/prepare_human_test_data.sh
   nextflow run main.nf -profile test,docker --reads 'test_data_human/fastq/*_R{1,2}.fastq.gz'
   ```

## Pipeline Components

### Quality Control

- **FastQC**: Quality assessment of raw and trimmed reads
- **MultiQC**: Aggregation of quality metrics across samples

### Read Trimming

- **Trimmomatic**: Removal of adapters and low-quality bases

### Alignment

- **STAR**: Splice-aware alignment for eukaryotic genomes
- **HISAT2**: Alternative splice-aware aligner

### Quantification

- **featureCounts**: Count reads mapped to genomic features
- **Salmon**: Transcript-level quantification with bias correction

### Differential Expression

- **DESeq2**: Statistical analysis of differential expression
- **Visualization**: MA plots, volcano plots, heatmaps, and PCA

### Reporting

- **R Markdown**: Generation of comprehensive HTML reports
- **Data Visualization**: Interactive plots and tables

## Results

What could you expect from using that pipeline:
   results -> fastQC/multiQC: 
      - this is quality control over the transcript and sequencing(*fastqc.html)
      - the adapters used in the sequencing(which you need to provide for yourself)(illumina_adapters.fa/illumina_adapters.txt)
   results -> featurecounts:
      - report about reads categorized in your BAM file(feature assignment)
      - gene lengths in the merge gene counts, only from one sample but could be set for more(reference from GTF/GFF annotation)
      - place in which those transcripts were located(on the genome - the chromosome, the start and the end coordintates, the strand orientation, length in bp's), the whole strand
      - summary of transcirpts in the table
   results -> star:
      - Aligned read files BAM file(s): One or more BAM files (e.g., Aligned.out.bam) containing the reads mapped to the reference genome.
      - Alignment summary logs:
         Log.final.out: A summary log detailing the alignment performance (e.g., percentage of uniquely mapped reads, number of reads processed).
         Log.out: A detailed log capturing runtime information, warnings, and additional performance metrics.
      - Splice junction information SJ.out.tab: A file listing all the splice junctions detected during the alignment process, which can be useful for downstream transcript isoform discovery or validation.
   results -> trimmed
      - illumina truqes adapter sequences files used for adapter trimming in RNA-seq

## Customization

### Configuration Files

The pipeline behavior can be customized by modifying the following files:

- `nextflow.config`: Main configuration file
- `conf/base.config`: Base configuration for resources
- `conf/docker.config`: Docker container configuration
- `conf/singularity.config`: Singularity container configuration

### Module Customization

Each module can be customized by modifying its main script:

- `modules/fastqc/main.nf`: FastQC module
- `modules/trimming/main.nf`: Trimming module
- `modules/alignment/star.nf`: STAR alignment module
- `modules/alignment/hisat2.nf`: HISAT2 alignment module
- `modules/quantification/featurecounts.nf`: featureCounts module
- `modules/quantification/salmon.nf`: Salmon module
- `modules/deseq2/main.nf`: DESeq2 module

## Troubleshooting

### Common Issues

1. **Insufficient resources**: Increase `--max_cpus` and `--max_memory` parameters
2. **Docker permission issues**: Ensure your user is in the docker group
3. **File not found errors**: Check that input paths are correct and files exist
4. **Nextflow version errors**: Update Nextflow to the required version

### Error Messages

| Error | Solution |
|-------|----------|
| `No such file or directory` | Check file paths and ensure files exist |
| `Docker pull failed` | Check internet connection and Docker installation |
| `Out of memory error` | Increase `--max_memory` parameter |
| `CPU usage 100%` | Increase `--max_cpus` parameter |

## Development

### Testing

The pipeline includes several test scripts:

```bash
# Run minimal test
bash bin/test_pipeline_mini.sh

# Run E. coli test
bash bin/test_pipeline_ecoli.sh

# Run synthetic data test
bash bin/test_pipeline_synthetic.sh

# Run human data test
bash bin/test_pipeline_human.sh
```

### Docker Containers

To build the Docker containers:

```bash
cd docker
bash build_containers.sh
```

## Citation
If you use this pipeline in your research, please cite:

```
Konrad Gorzelanczyk et al. (2025). A reproducible RNA-Seq analysis pipeline implemented in Nextflow.
GitHub: https://github.com/kGorze/RNA-Seq-Nextflow
```