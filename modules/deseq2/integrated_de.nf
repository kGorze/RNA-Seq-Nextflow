#!/usr/bin/env nextflow

/*
 * Integrated differential expression analysis workflow for RNA-Seq pipeline
 * Combines DESeq2 analysis with additional visualization and functional enrichment
 */

// Import DESeq2 module
include { DESEQ2_WORKFLOW } from '../deseq2/main'

workflow DE_ANALYSIS_WORKFLOW {
    take:
    counts_ch    // Channel with count matrices
    design_file  // Experimental design file
    contrasts_file // Contrasts definition file
    gtf_file     // Gene annotation GTF for gene name mapping
    
    main:
    // Run DESeq2 workflow
    DESEQ2_WORKFLOW(counts_ch, design_file, contrasts_file)
    
    // Generate sample correlation heatmap
    SAMPLE_CORRELATION(counts_ch, design_file)
    
    // Generate gene name mapping from GTF
    GENE_ID_MAPPING(gtf_file, DESEQ2_WORKFLOW.out.results)
    
    // Create integrated report
    INTEGRATED_DE_REPORT(
        DESEQ2_WORKFLOW.out.results,
        DESEQ2_WORKFLOW.out.plots,
        DESEQ2_WORKFLOW.out.report,
        SAMPLE_CORRELATION.out.correlation,
        GENE_ID_MAPPING.out.gene_mapping,
        design_file,
        contrasts_file
    )
    
    emit:
    results = DESEQ2_WORKFLOW.out.results
    plots = DESEQ2_WORKFLOW.out.plots
    report = INTEGRATED_DE_REPORT.out.report
    rdata = DESEQ2_WORKFLOW.out.rdata
}

process SAMPLE_CORRELATION {
    label 'process_low'
    
    container 'quay.io/biocontainers/r-base:4.1.0'
    
    publishDir "${params.outdir}/deseq2/correlation", mode: 'copy'
    
    input:
    path counts
    path design
    
    output:
    path "sample_correlation.pdf", emit: correlation
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    library(pheatmap)
    library(RColorBrewer)
    
    # Read count data
    counts_data <- read.csv("${counts}", row.names=1)
    
    # Remove length column if present
    if ("length" %in% colnames(counts_data)) {
        counts_data <- counts_data[, !colnames(counts_data) %in% "length"]
    }
    
    # Read experimental design
    design_data <- read.csv("${design}", row.names=1)
    
    # Ensure sample order in design matches count matrix
    design_data <- design_data[colnames(counts_data), , drop=FALSE]
    
    # Calculate correlation between samples
    sample_cor <- cor(counts_data, method="spearman")
    
    # Create annotation data frame
    anno <- data.frame(Condition = design_data\$condition)
    rownames(anno) <- colnames(counts_data)
    
    # Create correlation heatmap
    pdf("sample_correlation.pdf", width=10, height=8)
    pheatmap(sample_cor,
           main="Sample Correlation Heatmap",
           annotation_col=anno,
           annotation_row=anno,
           col=colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(255))
    dev.off()
    """
}

process GENE_ID_MAPPING {
    label 'process_low'
    
    container 'quay.io/biocontainers/bioconductor-rtracklayer:1.52.0--r41h399db7b_0'
    
    publishDir "${params.outdir}/deseq2/gene_mapping", mode: 'copy'
    
    input:
    path gtf
    path results
    
    output:
    path "gene_id_mapping.csv", emit: gene_mapping
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    library(rtracklayer)
    
    # Import GTF file
    gtf_data <- rtracklayer::import("${gtf}")
    
    # Extract gene information
    genes <- gtf_data[gtf_data\$type == "gene"]
    
    # Create mapping data frame
    gene_mapping <- data.frame(
        gene_id = genes\$gene_id,
        gene_name = ifelse(is.null(genes\$gene_name), genes\$gene_id, genes\$gene_name),
        gene_biotype = ifelse(is.null(genes\$gene_biotype), "unknown", genes\$gene_biotype)
    )
    
    # Remove duplicates
    gene_mapping <- gene_mapping[!duplicated(gene_mapping\$gene_id), ]
    
    # Write mapping to file
    write.csv(gene_mapping, file="gene_id_mapping.csv", row.names=FALSE)
    """
}

process INTEGRATED_DE_REPORT {
    label 'process_medium'
    
    container 'quay.io/biocontainers/bioconductor-deseq2:1.32.0--r41h399db7b_0'
    
    publishDir "${params.outdir}/deseq2/integrated_report", mode: 'copy'
    
    input:
    path results
    path plots
    path deseq2_report
    path correlation
    path gene_mapping
    path design
    path contrasts
    
    output:
    path "integrated_de_report.html", emit: report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Create integrated report template
    cat > integrated_de_report.Rmd << 'EOL'
---
title: "Integrated RNA-Seq Differential Expression Analysis Report"
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
```

# Overview

This report integrates all differential expression analysis results from the RNA-Seq pipeline.

## Experimental Design

```{r design}
# Read design file
design <- read.csv("${design}", row.names=1)

# Display as interactive table
DT::datatable(design, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Experimental Design")
```

## Contrasts Analyzed

```{r contrasts}
# Read contrasts file
contrasts <- read.csv("${contrasts}")

# Display as interactive table
DT::datatable(contrasts, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Contrasts Analyzed")
```

# Sample Correlation

The following heatmap shows the correlation between samples based on their gene expression profiles.

![Sample Correlation Heatmap](${correlation})

# Differential Expression Results

```{r read-results}
# Function to read and summarize results
summarize_results <- function(contrast_name) {
  # Read results file
  res_file <- paste0(contrast_name, "_results.csv")
  if (!file.exists(res_file)) {
    return(NULL)
  }
  
  res <- read.csv(res_file, row.names=1)
  
  # Count significant genes
  sig_up <- sum(res\$padj < 0.05 & res\$log2FoldChange > 0, na.rm=TRUE)
  sig_down <- sum(res\$padj < 0.05 & res\$log2FoldChange < 0, na.rm=TRUE)
  
  return(list(
    contrast = contrast_name,
    total = nrow(res),
    significant = sig_up + sig_down,
    up = sig_up,
    down = sig_down,
    results = res
  ))
}

# Read gene mapping
gene_mapping <- read.csv("${gene_mapping}")

# Process all contrasts
result_summaries <- list()
for (i in 1:nrow(contrasts)) {
  contrast_name <- contrasts[i, "name"]
  result_summaries[[contrast_name]] <- summarize_results(contrast_name)
}
```

## Summary of Differential Expression

The following table summarizes the number of differentially expressed genes for each contrast.

```{r summary-table}
# Create summary table
summary_df <- data.frame(
  Contrast = sapply(result_summaries, function(x) x\$contrast),
  Total_Genes = sapply(result_summaries, function(x) x\$total),
  Significant_Genes = sapply(result_summaries, function(x) x\$significant),
  Upregulated = sapply(result_summaries, function(x) x\$up),
  Downregulated = sapply(result_summaries, function(x) x\$down)
)

# Display as interactive table
DT::datatable(summary_df, 
              options = list(pageLength = 10),
              caption = "Summary of Differential Expression Results")
```

## Visualization of DE Results

```{r de-viz, fig.width=10, fig.height=6}
# Create bar plot of DE genes
ggplot(summary_df, aes(x = Contrast)) +
  geom_bar(aes(y = Upregulated, fill = "Upregulated"), stat = "identity") +
  geom_bar(aes(y = -Downregulated, fill = "Downregulated"), stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "darkred", "Downregulated" = "darkblue")) +
  labs(title = "Differentially Expressed Genes by Contrast",
       y = "Number of DE Genes",
       fill = "Regulation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

## Top Differentially Expressed Genes

```{r top-genes, results='asis'}
# Function to create interactive table of top DE genes with gene names
display_top_genes <- function(results, n=20) {
  # Sort by adjusted p-value
  results <- results[order(results\$padj), ]
  
  # Select top genes
  top_genes <- head(results, n)
  
  # Add gene names
  top_genes\$gene_id <- rownames(top_genes)
  top_genes <- merge(top_genes, gene_mapping, by="gene_id", all.x=TRUE)
  
  # Format table
  top_genes\$padj <- format(top_genes\$padj, digits=3, scientific=TRUE)
  top_genes\$pvalue <- format(top_genes\$pvalue, digits=3, scientific=TRUE)
  top_genes\$log2FoldChange <- round(top_genes\$log2FoldChange, 3)
  
  # Reorder columns
  top_genes <- top_genes[, c("gene_id", "gene_name", "gene_biotype", "log2FoldChange", "padj", "pvalue", "baseMean")]
  
  # Display as interactive table
  DT::datatable(top_genes, 
                options = list(pageLength = 10, scrollX = TRUE),
                caption = paste("Top", n, "Differentially Expressed Genes"))
}

# Display top genes for each contrast
for (i in 1:nrow(contrasts)) {
  contrast_name <- contrasts[i, "name"]
  control <- contrasts[i, "control"]
  treatment <- contrasts[i, "treatment"]
  
  cat("### ", contrast_name, " (", treatment, " vs ", control, ")\n\n", sep="")
  
  # Get results for this contrast
  results <- result_summaries[[contrast_name]]\$results
  
  # Display table
  display_top_genes(results)
  
  cat("\n\n")
}
```

# Detailed DESeq2 Report

For a more detailed analysis of the differential expression results, please refer to the [DESeq2 Report](${deseq2_report}).

# Conclusions

This report provides a comprehensive analysis of differential gene expression in the RNA-Seq dataset. The key findings are:

1. A total of `r sum(sapply(result_summaries, function(x) x\$significant))` genes were found to be differentially expressed across all contrasts.
2. The most significant contrast was `r summary_df\$Contrast[which.max(summary_df\$Significant_Genes)]` with `r max(summary_df\$Significant_Genes)` differentially expressed genes.
3. The sample correlation heatmap shows clear patterns of similarity between samples from the same condition.

## Next Steps

Based on these results, the following next steps are recommended:

1. Functional enrichment analysis of differentially expressed genes
2. Pathway analysis to identify affected biological processes
3. Validation of key findings using qPCR or other methods
4. Integration with other omics data if available

EOL

    # Render the R Markdown document
    rmarkdown::render("integrated_de_report.Rmd", output_file = "integrated_de_report.html")
    """
}
