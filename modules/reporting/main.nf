#!/usr/bin/env nextflow

/*
 * Comprehensive reporting module for RNA-Seq pipeline
 * Integrates results from all pipeline components into a cohesive final report
 */

workflow FINAL_REPORT_WORKFLOW {
    take:
    fastqc_results    // Channel with FastQC results
    multiqc_report    // Channel with MultiQC report
    alignment_stats   // Channel with alignment statistics
    quant_results     // Channel with quantification results
    de_results        // Channel with differential expression results
    de_report         // Channel with differential expression report
    pipeline_info     // Channel with pipeline execution information
    
    main:
    // Generate final HTML report
    GENERATE_FINAL_REPORT(
        fastqc_results,
        multiqc_report,
        alignment_stats,
        quant_results,
        de_results,
        de_report,
        pipeline_info
    )
    
    // Generate summary statistics
    GENERATE_SUMMARY_STATS(
        alignment_stats,
        quant_results,
        de_results
    )
    
    // Create results archive
    CREATE_RESULTS_ARCHIVE(
        GENERATE_FINAL_REPORT.out.report,
        GENERATE_SUMMARY_STATS.out.summary,
        fastqc_results,
        multiqc_report,
        alignment_stats,
        quant_results,
        de_results,
        de_report
    )
    
    emit:
    report = GENERATE_FINAL_REPORT.out.report
    summary = GENERATE_SUMMARY_STATS.out.summary
    archive = CREATE_RESULTS_ARCHIVE.out.archive
}

process GENERATE_FINAL_REPORT {
    label 'process_medium'
    
    container 'quay.io/biocontainers/r-rmarkdown:1.15--r40h6115d3f_1'
    
    publishDir "${params.outdir}/final_report", mode: 'copy'
    
    input:
    path fastqc_results
    path multiqc_report
    path alignment_stats
    path quant_results
    path de_results
    path de_report
    path pipeline_info
    
    output:
    path "final_report.html", emit: report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Create R Markdown template for final report
    cat > final_report.Rmd << 'EOL'
---
title: "RNA-Seq Analysis Final Report"
author: "Nextflow RNA-Seq Pipeline"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    highlight: tango
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(DT)
library(dplyr)
library(ggplot2)
library(knitr)
```

# RNA-Seq Analysis Overview

This report summarizes the results of the RNA-Seq analysis pipeline. The analysis was performed using a Nextflow pipeline that includes quality control, read alignment, transcript quantification, and differential expression analysis.

## Pipeline Information

```{r pipeline-info}
# Read pipeline info
pipeline_info <- readLines("${pipeline_info}")
```

- **Pipeline Version**: `r gsub("Pipeline version: ", "", grep("Pipeline version", pipeline_info, value=TRUE))`
- **Run Date**: `r gsub("Run date: ", "", grep("Run date", pipeline_info, value=TRUE))`
- **Execution Time**: `r gsub("Execution time: ", "", grep("Execution time", pipeline_info, value=TRUE))`

## Analysis Workflow

The RNA-Seq analysis pipeline consists of the following steps:

1. **Quality Control**: FastQC and MultiQC for assessing read quality
2. **Read Trimming**: Removal of adapters and low-quality bases
3. **Alignment**: Mapping reads to the reference genome
4. **Quantification**: Counting reads mapped to genes/transcripts
5. **Differential Expression**: Statistical analysis to identify differentially expressed genes

# Quality Control Results

The quality control step assessed the quality of the raw sequencing reads.

## MultiQC Report

The MultiQC report provides an aggregated view of quality control metrics across all samples.

<iframe src="${multiqc_report}" width="100%" height="600px"></iframe>

[Open MultiQC Report in New Tab](${multiqc_report})

# Alignment Results

```{r alignment-stats}
# Read alignment statistics
alignment_stats <- read.csv("${alignment_stats}", header=TRUE)

# Display as interactive table
DT::datatable(alignment_stats, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Alignment Statistics")
```

## Alignment Summary

```{r alignment-summary, fig.width=10, fig.height=6}
# Create bar plot of alignment rates
ggplot(alignment_stats, aes(x = Sample, y = Alignment_Rate, fill = Sample)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Alignment Rates by Sample",
       y = "Alignment Rate (%)",
       x = "Sample") +
  ylim(0, 100)
```

# Quantification Results

```{r quant-results}
# Read quantification summary
quant_summary <- read.csv("${quant_results}", header=TRUE)

# Display as interactive table
DT::datatable(quant_summary, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Quantification Summary")
```

## Gene Detection

```{r gene-detection, fig.width=10, fig.height=6}
# Create bar plot of detected genes
ggplot(quant_summary, aes(x = Sample, y = Detected_Genes, fill = Sample)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Detected Genes by Sample",
       y = "Number of Detected Genes",
       x = "Sample")
```

# Differential Expression Results

The differential expression analysis identified genes that are significantly differentially expressed between conditions.

## Differential Expression Summary

```{r de-summary}
# Read DE summary
de_summary <- read.csv("${de_results}", header=TRUE)

# Display as interactive table
DT::datatable(de_summary, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Differential Expression Summary")
```

## Visualization of DE Results

```{r de-viz, fig.width=10, fig.height=6}
# Create bar plot of DE genes
ggplot(de_summary, aes(x = Contrast)) +
  geom_bar(aes(y = Upregulated, fill = "Upregulated"), stat = "identity") +
  geom_bar(aes(y = -Downregulated, fill = "Downregulated"), stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "darkred", "Downregulated" = "darkblue")) +
  labs(title = "Differentially Expressed Genes by Contrast",
       y = "Number of DE Genes",
       fill = "Regulation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

# Detailed Reports

For more detailed results, please refer to the following reports:

- [MultiQC Report](${multiqc_report})
- [Differential Expression Report](${de_report})

# Conclusions

This report provides an overview of the RNA-Seq analysis results. The key findings are:

1. The quality control metrics indicate that the sequencing data is of good quality.
2. The alignment rates are high, suggesting successful mapping to the reference genome.
3. A substantial number of genes were detected across all samples.
4. The differential expression analysis identified significant changes in gene expression between conditions.

## Next Steps

Based on these results, the following next steps are recommended:

1. Functional enrichment analysis of differentially expressed genes
2. Pathway analysis to identify affected biological processes
3. Validation of key findings using qPCR or other methods
4. Integration with other omics data if available

EOL

    # Render the R Markdown document
    rmarkdown::render("final_report.Rmd", output_file = "final_report.html")
    """
}

process GENERATE_SUMMARY_STATS {
    label 'process_low'
    
    container 'quay.io/biocontainers/r-base:4.1.0'
    
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path alignment_stats
    path quant_results
    path de_results
    
    output:
    path "summary_statistics.csv", emit: summary
    path "summary_statistics.json", emit: json
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Read input files
    alignment_stats <- read.csv("${alignment_stats}", header=TRUE)
    quant_results <- read.csv("${quant_results}", header=TRUE)
    de_results <- read.csv("${de_results}", header=TRUE)
    
    # Calculate summary statistics
    summary_stats <- data.frame(
        Metric = c(
            "Total Samples",
            "Mean Alignment Rate (%)",
            "Mean Detected Genes",
            "Total Contrasts Analyzed",
            "Total DE Genes",
            "Mean DE Genes per Contrast"
        ),
        Value = c(
            nrow(alignment_stats),
            mean(alignment_stats\$Alignment_Rate),
            mean(quant_results\$Detected_Genes),
            nrow(de_results),
            sum(de_results\$Significant_Genes),
            mean(de_results\$Significant_Genes)
        )
    )
    
    # Write summary statistics to CSV
    write.csv(summary_stats, file="summary_statistics.csv", row.names=FALSE)
    
    # Write summary statistics to JSON
    json_data <- paste0(
        "{\n",
        paste0("  \"", summary_stats\$Metric, "\": ", summary_stats\$Value, collapse=",\n"),
        "\n}"
    )
    writeLines(json_data, "summary_statistics.json")
    """
}

process CREATE_RESULTS_ARCHIVE {
    label 'process_low'
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path final_report
    path summary_stats
    path fastqc_results
    path multiqc_report
    path alignment_stats
    path quant_results
    path de_results
    path de_report
    
    output:
    path "results_archive.zip", emit: archive
    
    script:
    """
    # Create directory structure
    mkdir -p results/qc
    mkdir -p results/alignment
    mkdir -p results/quantification
    mkdir -p results/de_analysis
    mkdir -p results/reports
    
    # Copy files to appropriate directories
    cp ${final_report} results/reports/
    cp ${summary_stats} results/
    cp ${fastqc_results} results/qc/
    cp ${multiqc_report} results/qc/
    cp ${alignment_stats} results/alignment/
    cp ${quant_results} results/quantification/
    cp ${de_results} results/de_analysis/
    cp ${de_report} results/de_analysis/
    
    # Create README file
    cat > results/README.txt << EOL
RNA-Seq Analysis Results
========================

This archive contains the results of the RNA-Seq analysis pipeline.

Directory Structure:
-------------------
- qc/: Quality control results
- alignment/: Alignment statistics
- quantification/: Gene/transcript quantification results
- de_analysis/: Differential expression analysis results
- reports/: Summary reports

Main Files:
----------
- reports/final_report.html: Comprehensive analysis report
- summary_statistics.csv: Summary statistics of the analysis

For more information, please refer to the documentation.
EOL
    
    # Create ZIP archive
    zip -r results_archive.zip results/
    """
}
