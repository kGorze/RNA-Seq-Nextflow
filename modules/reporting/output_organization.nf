#!/usr/bin/env nextflow

/*
 * Output organization module for RNA-Seq pipeline
 * Organizes and structures output files for better usability
 */

workflow OUTPUT_ORGANIZATION_WORKFLOW {
    take:
    fastqc_results    // Channel with FastQC results
    multiqc_results   // Channel with MultiQC results
    trimming_results  // Channel with trimming results
    alignment_results // Channel with alignment results
    quant_results     // Channel with quantification results
    de_results        // Channel with differential expression results
    reports           // Channel with reports
    
    main:
    // Organize output files
    ORGANIZE_OUTPUTS(
        fastqc_results,
        multiqc_results,
        trimming_results,
        alignment_results,
        quant_results,
        de_results,
        reports
    )
    
    // Generate output documentation
    GENERATE_OUTPUT_DOCUMENTATION()
    
    emit:
    documentation = GENERATE_OUTPUT_DOCUMENTATION.out.documentation
}

process ORGANIZE_OUTPUTS {
    label 'process_low'
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path fastqc_results
    path multiqc_results
    path trimming_results
    path alignment_results
    path quant_results
    path de_results
    path reports
    
    output:
    path "organized_outputs", emit: organized_outputs
    
    script:
    """
    # Create organized directory structure
    mkdir -p organized_outputs/01_quality_control
    mkdir -p organized_outputs/02_trimming
    mkdir -p organized_outputs/03_alignment
    mkdir -p organized_outputs/04_quantification
    mkdir -p organized_outputs/05_differential_expression
    mkdir -p organized_outputs/06_reports
    
    # Copy files to appropriate directories
    cp -r ${fastqc_results} organized_outputs/01_quality_control/
    cp -r ${multiqc_results} organized_outputs/01_quality_control/
    cp -r ${trimming_results} organized_outputs/02_trimming/
    cp -r ${alignment_results} organized_outputs/03_alignment/
    cp -r ${quant_results} organized_outputs/04_quantification/
    cp -r ${de_results} organized_outputs/05_differential_expression/
    cp -r ${reports} organized_outputs/06_reports/
    
    # Create README file
    cat > organized_outputs/README.md << EOL
# RNA-Seq Pipeline Output Directory

This directory contains the organized outputs from the RNA-Seq pipeline.

## Directory Structure

- **01_quality_control/**: FastQC and MultiQC results for raw and trimmed reads
- **02_trimming/**: Trimmed reads and trimming logs
- **03_alignment/**: Alignment files (BAM) and alignment statistics
- **04_quantification/**: Gene and transcript quantification results
- **05_differential_expression/**: Differential expression analysis results
- **06_reports/**: Summary reports and visualizations

## Key Files

- **06_reports/final_report.html**: Comprehensive analysis report
- **06_reports/summary_statistics.csv**: Summary statistics of the analysis
- **01_quality_control/multiqc/multiqc_report.html**: Quality control summary
- **05_differential_expression/deseq2/integrated_de_report.html**: Detailed differential expression results

For more information, please refer to the pipeline documentation.
EOL
    """
}

process GENERATE_OUTPUT_DOCUMENTATION {
    label 'process_low'
    
    publishDir "${params.outdir}/documentation", mode: 'copy'
    
    output:
    path "output_documentation.html", emit: documentation
    
    script:
    """
    # Generate HTML documentation
    cat > output_documentation.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RNA-Seq Pipeline Output Documentation</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }
        h1, h2, h3 {
            color: #2c3e50;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        table {
            border-collapse: collapse;
            width: 100%;
            margin-bottom: 20px;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        th {
            background-color: #f2f2f2;
        }
        tr:nth-child(even) {
            background-color: #f9f9f9;
        }
        .file-path {
            font-family: monospace;
            background-color: #f5f5f5;
            padding: 2px 4px;
            border-radius: 3px;
        }
        .note {
            background-color: #e7f3fe;
            border-left: 6px solid #2196F3;
            padding: 10px;
            margin-bottom: 15px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>RNA-Seq Pipeline Output Documentation</h1>
        
        <div class="note">
            <p><strong>Note:</strong> This document describes the output files generated by the RNA-Seq pipeline and provides guidance on interpreting the results.</p>
        </div>
        
        <h2>Directory Structure</h2>
        
        <p>The pipeline organizes outputs into the following directory structure:</p>
        
        <table>
            <tr>
                <th>Directory</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><span class="file-path">01_quality_control/</span></td>
                <td>FastQC and MultiQC results for raw and trimmed reads</td>
            </tr>
            <tr>
                <td><span class="file-path">02_trimming/</span></td>
                <td>Trimmed reads and trimming logs</td>
            </tr>
            <tr>
                <td><span class="file-path">03_alignment/</span></td>
                <td>Alignment files (BAM) and alignment statistics</td>
            </tr>
            <tr>
                <td><span class="file-path">04_quantification/</span></td>
                <td>Gene and transcript quantification results</td>
            </tr>
            <tr>
                <td><span class="file-path">05_differential_expression/</span></td>
                <td>Differential expression analysis results</td>
            </tr>
            <tr>
                <td><span class="file-path">06_reports/</span></td>
                <td>Summary reports and visualizations</td>
            </tr>
        </table>
        
        <h2>Key Output Files</h2>
        
        <h3>Quality Control</h3>
        
        <table>
            <tr>
                <th>File</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><span class="file-path">01_quality_control/fastqc/*_fastqc.html</span></td>
                <td>FastQC reports for individual samples</td>
            </tr>
            <tr>
                <td><span class="file-path">01_quality_control/multiqc/multiqc_report.html</span></td>
                <td>Aggregated quality control metrics across all samples</td>
            </tr>
        </table>
        
        <h3>Trimming</h3>
        
        <table>
            <tr>
                <th>File</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><span class="file-path">02_trimming/*_trimmed.fastq.gz</span></td>
                <td>Trimmed read files</td>
            </tr>
            <tr>
                <td><span class="file-path">02_trimming/*_trimming_report.txt</span></td>
                <td>Trimming statistics and log files</td>
            </tr>
        </table>
        
        <h3>Alignment</h3>
        
        <table>
            <tr>
                <th>File</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><span class="file-path">03_alignment/star/*.sorted.bam</span></td>
                <td>Sorted BAM files from STAR alignment</td>
            </tr>
            <tr>
                <td><span class="file-path">03_alignment/hisat2/*.sorted.bam</span></td>
                <td>Sorted BAM files from HISAT2 alignment</td>
            </tr>
            <tr>
                <td><span class="file-path">03_alignment/*/stats/*.stats.txt</span></td>
                <td>Alignment statistics</td>
            </tr>
        </table>
        
        <h3>Quantification</h3>
        
        <table>
            <tr>
                <th>File</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><span class="file-path">04_quantification/featurecounts/merged/merged_gene_counts.csv</span></td>
                <td>Merged gene-level count matrix from featureCounts</td>
            </tr>
            <tr>
                <td><span class="file-path">04_quantification/salmon/tximport/gene_counts.csv</span></td>
                <td>Gene-level count matrix from Salmon</td>
            </tr>
            <tr>
                <td><span class="file-path">04_quantification/salmon/quant/*/quant.sf</span></td>
                <td>Transcript-level quantification from Salmon</td>
            </tr>
        </table>
        
        <h3>Differential Expression</h3>
        
        <table>
            <tr>
                <th>File</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><span class="file-path">05_differential_expression/deseq2/results/*_results.csv</span></td>
                <td>Differential expression results for each contrast</td>
            </tr>
            <tr>
                <td><span class="file-path">05_differential_expression/deseq2/results/*.pdf</span></td>
                <td>Visualizations of differential expression results (MA plots, volcano plots, heatmaps)</td>
            </tr>
            <tr>
                <td><span class="file-path">05_differential_expression/deseq2/report/deseq2_report.html</span></td>
                <td>Detailed DESeq2 analysis report</td>
            </tr>
            <tr>
                <td><span class="file-path">05_differential_expression/deseq2/integrated_report/integrated_de_report.html</span></td>
                <td>Integrated differential expression analysis report</td>
            </tr>
        </table>
        
        <h3>Reports</h3>
        
        <table>
            <tr>
                <th>File</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><span class="file-path">06_reports/final_report.html</span></td>
                <td>Comprehensive analysis report integrating all pipeline components</td>
            </tr>
            <tr>
                <td><span class="file-path">06_reports/summary_statistics.csv</span></td>
                <td>Summary statistics of the analysis</td>
            </tr>
            <tr>
                <td><span class="file-path">06_reports/pipeline_info.txt</span></td>
                <td>Pipeline execution information</td>
            </tr>
        </table>
        
        <h2>Interpreting Results</h2>
        
        <h3>Quality Control</h3>
        
        <p>The FastQC and MultiQC reports provide information about the quality of the sequencing data. Key metrics to look for include:</p>
        
        <ul>
            <li>Per base sequence quality (should be high, typically >30)</li>
            <li>Adapter content (should be low or absent after trimming)</li>
            <li>Sequence duplication levels (some duplication is expected in RNA-Seq)</li>
            <li>GC content distribution (should follow a normal distribution)</li>
        </ul>
        
        <h3>Alignment</h3>
        
        <p>The alignment statistics provide information about how well the reads mapped to the reference genome. Key metrics include:</p>
        
        <ul>
            <li>Alignment rate (typically >80% for good quality RNA-Seq data)</li>
            <li>Uniquely mapped reads (higher is better)</li>
            <li>Multi-mapped reads (some multi-mapping is expected for RNA-Seq)</li>
        </ul>
        
        <h3>Differential Expression</h3>
        
        <p>The differential expression results identify genes that are significantly differentially expressed between conditions. Key columns in the results files include:</p>
        
        <ul>
            <li><strong>log2FoldChange</strong>: Log2 fold change in expression between conditions</li>
            <li><strong>padj</strong>: Adjusted p-value (FDR) for statistical significance</li>
            <li><strong>baseMean</strong>: Average expression level across all samples</li>
        </ul>
        
        <p>Genes with <strong>padj < 0.05</strong> and <strong>|log2FoldChange| > 1</strong> are typically considered significantly differentially expressed.</p>
        
        <h2>Next Steps</h2>
        
        <p>After reviewing the results, you may want to:</p>
        
        <ol>
            <li>Perform functional enrichment analysis of differentially expressed genes</li>
            <li>Conduct pathway analysis to identify affected biological processes</li>
            <li>Validate key findings using qPCR or other methods</li>
            <li>Integrate with other omics data if available</li>
        </ol>
        
        <div class="note">
            <p><strong>Note:</strong> For any questions or issues with the pipeline outputs, please refer to the pipeline documentation or contact the pipeline developers.</p>
        </div>
    </div>
</body>
</html>
EOL
    """
}
