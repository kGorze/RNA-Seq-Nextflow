#!/usr/bin/env nextflow

/*
 * Enhanced DESeq2 differential expression analysis workflow for RNA-Seq pipeline
 * Includes comprehensive statistical analysis and visualization
 */

workflow DESEQ2_WORKFLOW {
    take:
    counts_ch    // Channel with count matrices
    design_file  // Experimental design file
    contrasts_file // Contrasts definition file
    
    main:
    // Validate design and contrasts files
    VALIDATE_DESIGN_CONTRASTS(design_file, contrasts_file)
    
    // Run DESeq2 analysis
    DESEQ2_ANALYSIS(counts_ch, design_file, contrasts_file)
    
    // Generate comprehensive report
    DESEQ2_REPORT(
        DESEQ2_ANALYSIS.out.results,
        DESEQ2_ANALYSIS.out.normalized_counts,
        DESEQ2_ANALYSIS.out.rdata,
        design_file,
        contrasts_file
    )
    
    emit:
    results = DESEQ2_ANALYSIS.out.results
    plots = DESEQ2_ANALYSIS.out.plots
    report = DESEQ2_REPORT.out.report
    rdata = DESEQ2_ANALYSIS.out.rdata
}

process VALIDATE_DESIGN_CONTRASTS {
    label 'process_low'
    
    container 'quay.io/biocontainers/r-base:4.1.0'
    
    input:
    path design
    path contrasts
    
    output:
    path "design_validation.txt", emit: validation
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Read design file
    design_data <- read.csv("${design}", row.names=1)
    
    # Read contrasts file
    contrasts_data <- read.csv("${contrasts}")
    
    # Validate design file
    sink("design_validation.txt")
    cat("Design file validation\\n")
    cat("---------------------\\n")
    
    # Check if design file has required columns
    if (!"condition" %in% colnames(design_data)) {
        cat("ERROR: Design file must contain a 'condition' column\\n")
        quit(status=1)
    }
    
    cat("Design file contains", nrow(design_data), "samples\\n")
    cat("Conditions found:", paste(unique(design_data\$condition), collapse=", "), "\\n\\n")
    
    # Validate contrasts file
    cat("Contrasts file validation\\n")
    cat("------------------------\\n")
    
    # Check if contrasts file has required columns
    required_cols <- c("name", "control", "treatment")
    missing_cols <- required_cols[!required_cols %in% colnames(contrasts_data)]
    
    if (length(missing_cols) > 0) {
        cat("ERROR: Contrasts file missing required columns:", paste(missing_cols, collapse=", "), "\\n")
        quit(status=1)
    }
    
    cat("Contrasts file contains", nrow(contrasts_data), "comparisons\\n")
    
    # Check if all conditions in contrasts file exist in design file
    all_conditions <- unique(c(contrasts_data\$control, contrasts_data\$treatment))
    missing_conditions <- all_conditions[!all_conditions %in% unique(design_data\$condition)]
    
    if (length(missing_conditions) > 0) {
        cat("ERROR: The following conditions in contrasts file are not found in design file:", 
            paste(missing_conditions, collapse=", "), "\\n")
        quit(status=1)
    }
    
    cat("All conditions in contrasts file exist in design file\\n")
    cat("Validation successful\\n")
    sink()
    """
}

process DESEQ2_ANALYSIS {
    label 'process_medium'
    
    container 'quay.io/biocontainers/bioconductor-deseq2:1.32.0--r41h399db7b_0'
    
    publishDir "${params.outdir}/deseq2/results", mode: 'copy'
    
    input:
    path counts
    path design
    path contrasts
    
    output:
    path "*.csv", emit: results
    path "*.pdf", emit: plots
    path "*.RData", emit: rdata
    path "normalized_counts.csv", emit: normalized_counts
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(EnhancedVolcano)
    
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
    
    # Convert condition to factor
    design_data\$condition <- factor(design_data\$condition)
    
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
        countData = round(as.matrix(counts_data)),
        colData = design_data,
        design = ~ condition
    )
    
    # Filter low count genes
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    
    # Run DESeq2 analysis
    dds <- DESeq(dds)
    
    # Save normalized counts
    normalized_counts <- counts(dds, normalized=TRUE)
    write.csv(normalized_counts, file="normalized_counts.csv")
    
    # Read contrasts
    contrasts_data <- read.csv("${contrasts}")
    
    # Generate results for each contrast
    for (i in 1:nrow(contrasts_data)) {
        contrast_name <- contrasts_data[i, "name"]
        control <- contrasts_data[i, "control"]
        treatment <- contrasts_data[i, "treatment"]
        
        # Get results
        res <- results(dds, contrast=c("condition", treatment, control))
        res <- res[order(res\$padj), ]
        
        # Write results to file
        write.csv(as.data.frame(res), file=paste0(contrast_name, "_results.csv"))
        
        # Create MA plot
        pdf(paste0(contrast_name, "_MA_plot.pdf"))
        plotMA(res, main=paste0(contrast_name, " MA Plot"))
        dev.off()
        
        # Create volcano plot
        pdf(paste0(contrast_name, "_volcano_plot.pdf"), width=10, height=8)
        EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = paste0(contrast_name, ' Volcano Plot'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      pointSize = 3.0,
                      labSize = 3.0)
        dev.off()
        
        # Create heatmap of top DE genes
        pdf(paste0(contrast_name, "_heatmap.pdf"), width=10, height=12)
        # Get top DE genes
        de_genes <- rownames(res)[res\$padj < 0.05 & abs(res\$log2FoldChange) > 1]
        if (length(de_genes) > 50) {
            de_genes <- de_genes[1:50]  # Limit to top 50 genes
        }
        
        if (length(de_genes) > 0) {
            # Extract normalized counts for DE genes
            de_counts <- normalized_counts[de_genes, ]
            
            # Scale counts for heatmap
            de_counts_scaled <- t(scale(t(de_counts)))
            
            # Create annotation data frame
            anno <- data.frame(Condition = design_data\$condition)
            rownames(anno) <- colnames(de_counts)
            
            # Create heatmap
            pheatmap(de_counts_scaled,
                   annotation_col = anno,
                   main = paste0(contrast_name, " - Top DE Genes"),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   fontsize_row = 8,
                   fontsize_col = 8)
        } else {
            plot(1, type="n", axes=FALSE, xlab="", ylab="")
            text(1, 1, "No significant DE genes found")
        }
        dev.off()
    }
    
    # Create PCA plot
    vsd <- vst(dds, blind=FALSE)
    pdf("PCA_plot.pdf", width=10, height=8)
    pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle("PCA Plot") +
        theme_minimal() +
        theme(legend.position="right")
    dev.off()
    
    # Create sample distance heatmap
    pdf("sample_distance_heatmap.pdf", width=10, height=8)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(vsd)
    colnames(sampleDistMatrix) <- colnames(vsd)
    pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           main="Sample Distance Heatmap",
           annotation_col=data.frame(Condition=design_data\$condition, row.names=rownames(design_data)),
           col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
    dev.off()
    
    # Save DESeq2 object for further analysis
    save(dds, vsd, file="deseq2_analysis.RData")
    """
}

process DESEQ2_REPORT {
    label 'process_medium'
    
    container 'quay.io/biocontainers/bioconductor-deseq2:1.32.0--r41h399db7b_0'
    
    publishDir "${params.outdir}/deseq2/report", mode: 'copy'
    
    input:
    path results
    path normalized_counts
    path rdata
    path design
    path contrasts
    
    output:
    path "deseq2_report.html", emit: report
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Create R Markdown template
    cat > deseq2_report.Rmd << 'EOL'
---
title: "RNA-Seq Differential Expression Analysis Report"
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
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(DT)
library(plotly)
```

# Overview

This report summarizes the differential expression analysis results from the RNA-Seq pipeline.

## Experimental Design

```{r design}
# Read design file
design <- read.csv("${design}", row.names=1)

# Display as interactive table
DT::datatable(design, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Experimental Design")
```

## Analysis Summary

```{r load-data}
# Load DESeq2 data
load("${rdata}")

# Read contrasts
contrasts <- read.csv("${contrasts}")
```

The differential expression analysis was performed using DESeq2 on `r nrow(design)` samples across `r length(unique(design\$condition))` conditions. A total of `r nrow(dds)` genes were analyzed after filtering low-count genes.

# Quality Control

## PCA Plot

The Principal Component Analysis (PCA) plot shows the overall similarity between samples based on their gene expression profiles.

```{r pca, fig.width=10, fig.height=8}
# Create PCA plot
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal() +
  theme(legend.position="right")

# Convert to plotly for interactivity
ggplotly(p)
```

## Sample Distance Heatmap

The sample distance heatmap shows the similarity between samples based on their gene expression profiles.

```{r sample-dist, fig.width=10, fig.height=8}
# Create sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

# Create annotation data frame
anno <- data.frame(Condition = design\$condition)
rownames(anno) <- rownames(design)

pheatmap(sampleDistMatrix,
       clustering_distance_rows=sampleDists,
       clustering_distance_cols=sampleDists,
       main="Sample Distance Heatmap",
       annotation_col=anno,
       col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
```

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

## Detailed Results by Contrast

```{r results-by-contrast, results='asis'}
# Function to create interactive table of top DE genes
display_top_genes <- function(results, n=20) {
  # Sort by adjusted p-value
  results <- results[order(results\$padj), ]
  
  # Select top genes
  top_genes <- head(results, n)
  
  # Format table
  top_genes\$padj <- format(top_genes\$padj, digits=3, scientific=TRUE)
  top_genes\$pvalue <- format(top_genes\$pvalue, digits=3, scientific=TRUE)
  top_genes\$log2FoldChange <- round(top_genes\$log2FoldChange, 3)
  
  # Display as interactive table
  DT::datatable(top_genes, 
                options = list(pageLength = 10, scrollX = TRUE),
                caption = paste("Top", n, "Differentially Expressed Genes"))
}

# Display results for each contrast
for (i in 1:nrow(contrasts)) {
  contrast_name <- contrasts[i, "name"]
  control <- contrasts[i, "control"]
  treatment <- contrasts[i, "treatment"]
  
  cat("### ", contrast_name, " (", treatment, " vs ", control, ")\n\n", sep="")
  
  # Display MA plot
  cat("#### MA Plot\n\n")
  cat("![", contrast_name, " MA Plot](", contrast_name, "_MA_plot.pdf)\n\n", sep="")
  
  # Display volcano plot
  cat("#### Volcano Plot\n\n")
  cat("![", contrast_name, " Volcano Plot](", contrast_name, "_volcano_plot.pdf)\n\n", sep="")
  
  # Display heatmap
  cat("#### Heatmap of Top DE Genes\n\n")
  cat("![", contrast_name, " Heatmap](", contrast_name, "_heatmap.pdf)\n\n", sep="")
  
  # Display top DE genes table
  cat("#### Top Differentially Expressed Genes\n\n")
  
  # Get results for this contrast
  results <- result_summaries[[contrast_name]]\$results
  
  # Display table
  display_top_genes(results)
  
  cat("\n\n")
}
```

# Gene Expression Patterns

## Normalized Counts

The following table shows the normalized counts for all samples.

```{r normalized-counts}
# Read normalized counts
norm_counts <- read.csv("${normalized_counts}", row.names=1)

# Display as interactive table (first 20 genes)
DT::datatable(head(norm_counts, 20), 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Normalized Counts (First 20 Genes)")
```

# Conclusions

This report provides a comprehensive analysis of differential gene expression in the RNA-Seq dataset. The key findings are:

1. The PCA plot shows clear separation between sample groups, indicating distinct gene expression patterns between conditions.
2. A total of `r sum(sapply(result_summaries, function(x) x\$significant))` genes were found to be differentially expressed across all contrasts.
3. The most significant contrast was `r summary_df\$Contrast[which.max(summary_df\$Significant_Genes)]` with `r max(summary_df\$Significant_Genes)` differentially expressed genes.

## Next Steps

Based on these results, the following next steps are recommended:

1. Functional enrichment analysis of differentially expressed genes
2. Pathway analysis to identify affected biological processes
3. Validation of key findings using qPCR or other methods

EOL

    # Render the R Markdown document
    rmarkdown::render("deseq2_report.Rmd", output_file = "deseq2_report.html")
    """
}
