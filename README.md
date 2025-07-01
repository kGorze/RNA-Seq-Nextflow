# RNA-Seq Nextflow Pipeline

A reproducible, containerized RNA‑seq workflow from FASTQ to differential expression implemented in Nextflow.
![image](https://github.com/user-attachments/assets/610a28da-9498-4bcb-aa79-8f7948121221)


## Quick start

```bash
git clone https://github.com/kGorze/RNA-Seq-Nextflow
cd RNA-Seq-Nextflow
nextflow run main.nf -profile test,docker

# run on your data
nextflow run main.nf -profile docker \
  --reads 'reads/*_R{1,2}.fastq.gz' \
  --genome genome.fa --gtf genes.gtf \
  --outdir results
```

## Minimal parameters

| flag       | meaning                                 |
| ---------- | --------------------------------------- |
| `--reads`  | Paired FASTQ glob (`*_R{1,2}.fastq.gz`) |
| `--genome` | Reference genome FASTA                  |
| `--gtf`    | Gene annotation GTF                     |
| `--design` | Sample design CSV for DESeq2 (optional) |

Full list: `docs/parameters.md`.

## Profiles

* **docker** / **singularity** – containerised (default docker)
* **test** – tiny dataset for CI / validation
* **slurm** / **sge** / **standard** – customise for your cluster / local install

Combine with commas: `-profile test,docker`.

## Outputs (key files)

| path                                       | content           |
| ------------------------------------------ | ----------------- |
| `01_quality_control/multiqc_report.html`   | QC summary        |
| `03_alignment/*/Log.final.out`             | Alignment stats   |
| `04_quantification/merged_gene_counts.csv` | Gene counts       |
| `05_differential_expression/deseq2/*.csv`  | DE results        |
| `06_reports/final_report.html`             | End‑to‑end report |


## Cite

> Gorzelanczyk K. (2025) Reproducible RNA‑Seq analysis with Nextflow. [https://github.com/kGorze/RNA-Seq-Nextflow](https://github.com/kGorze/RNA-Seq-Nextflow)
