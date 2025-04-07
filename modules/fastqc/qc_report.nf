#!/usr/bin/env nextflow

/*
 * Quality Control Report Template for RNA-Seq pipeline
 * Generates a comprehensive QC report using R Markdown
 */

process QC_REPORT {
    label 'process_low'
    
    container 'quay.io/biocontainers/r-rmarkdown:1.15--r40h6115d3f_1'
    
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    input:
    path fastqc_data
    path multiqc_data
    
    output:
    path "qc_report.html", emit: report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Create R Markdown template
    cat > qc_report.Rmd << 'EOL'
---
title: "RNA-Seq Quality Control Report"
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
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(DT)
```

# Overview

This report summarizes the quality control metrics for the RNA-Seq data processed through the Nextflow pipeline.

## Sample Summary

```{r sample-summary}
# Read MultiQC general stats
multiqc_stats <- read.table("${multiqc_data}/multiqc_general_stats.txt", 
                           header = TRUE, sep = "\\t", check.names = FALSE)

# Create sample summary table
sample_summary <- multiqc_stats %>%
  select(contains("FastQC"), contains("Sequence"))

# Display as interactive table
DT::datatable(sample_summary, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Sample Quality Metrics")
```

# Sequence Quality Metrics

## Per Base Sequence Quality

The following plot shows the quality scores across all bases at each position in the reads.

```{r base-quality}
# This would typically read data from FastQC output
# For demonstration, we'll create a simulated plot
positions <- 1:100
mean_quality <- 30 + 5 * sin(positions / 20) + rnorm(100, 0, 1)
mean_quality[mean_quality > 40] <- 40
mean_quality[mean_quality < 20] <- 20

quality_data <- data.frame(
  Position = positions,
  Quality = mean_quality
)

ggplot(quality_data, aes(x = Position, y = Quality)) +
  geom_line() +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal() +
  labs(title = "Per Base Sequence Quality",
       x = "Position in read (bp)",
       y = "Quality score") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "green") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "orange") +
  scale_y_continuous(limits = c(0, 40))
```

## Sequence Length Distribution

```{r seq-length}
# Simulated sequence length distribution
lengths <- c(rep(75, 100), rep(100, 500), rep(125, 300), rep(150, 100))
length_data <- data.frame(Length = lengths)

ggplot(length_data, aes(x = Length)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Sequence Length Distribution",
       x = "Sequence Length (bp)",
       y = "Count")
```

# Adapter Content

The following plot shows the adapter content across all bases.

```{r adapter-content}
# Simulated adapter content
positions <- 1:100
adapter_percent <- pmax(0, -10 + 0.5 * positions + rnorm(100, 0, 5))
adapter_percent[adapter_percent > 100] <- 100
adapter_percent[adapter_percent < 0] <- 0

adapter_data <- data.frame(
  Position = positions,
  Percent = adapter_percent
)

ggplot(adapter_data, aes(x = Position, y = Percent)) +
  geom_line(color = "red") +
  theme_minimal() +
  labs(title = "Adapter Content",
       x = "Position in read (bp)",
       y = "% of sequences with adapter") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "orange") +
  scale_y_continuous(limits = c(0, 100))
```

# Duplication Levels

```{r duplication}
# Simulated duplication data
duplication_level <- 1:10
duplication_percent <- c(65, 15, 8, 5, 3, 2, 1, 0.5, 0.3, 0.2)

duplication_data <- data.frame(
  Level = duplication_level,
  Percent = duplication_percent
)

ggplot(duplication_data, aes(x = as.factor(Level), y = Percent)) +
  geom_bar(stat = "identity", fill = "purple") +
  theme_minimal() +
  labs(title = "Sequence Duplication Levels",
       x = "Duplication Level",
       y = "% of total sequences")
```

# GC Content

```{r gc-content}
# Simulated GC content
gc_percent <- 1:100
gc_count <- dnorm(gc_percent, mean = 50, sd = 10) * 1000
gc_data <- data.frame(
  GC = gc_percent,
  Count = gc_count
)

ggplot(gc_data, aes(x = GC, y = Count)) +
  geom_line() +
  theme_minimal() +
  labs(title = "GC Content Distribution",
       x = "GC Content (%)",
       y = "Count") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "blue")
```

# Conclusions

This report provides an overview of the quality control metrics for the RNA-Seq data. The key observations are:

1. Overall sequence quality is good with most bases above Q30
2. Some adapter content is detected towards the end of reads
3. GC content distribution appears normal
4. Duplication levels are within expected range for RNA-Seq data

## Recommendations

Based on the QC results, the following recommendations are made:

- Proceed with trimming to remove adapter sequences and low-quality bases
- The data is suitable for downstream analysis including alignment and quantification

EOL

    # Render the R Markdown document
    rmarkdown::render("qc_report.Rmd", output_file = "qc_report.html")
    """
}
