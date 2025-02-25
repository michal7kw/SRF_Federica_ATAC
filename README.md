# ATAC-seq Analysis Pipeline

This Snakemake pipeline processes ATAC-seq data and performs differential analysis between experimental conditions.

## Prerequisites

The pipeline requires the following software to be installed in the conda environment named "snakemake":

- Snakemake (≥6.0)
- FastQC
- Trimmomatic
- Bowtie2
- Samtools
- Picard
- MACS2
- R with the following packages:
  - DiffBind
  - tidyverse

## Directory Structure

```
.
├── config.yaml           # Configuration file
├── Snakefile            # Pipeline rules
├── scripts/
│   └── run_diffbind.R   # R script for differential analysis
└── results/             # Output directory (created during execution)
    ├── qc/              # Quality control results
    ├── trimmed/         # Trimmed fastq files
    ├── bam/             # Aligned and processed BAM files
    ├── peaks/           # MACS2 peak calls
    └── diffbind/        # Differential binding results
```

## Configuration

1. Edit `config.yaml` to specify:
   - Sample information
   - Reference genome paths
   - Tool parameters
   - Differential analysis comparisons

2. Required reference files:
   - Genome FASTA file
   - Bowtie2 index
   - Blacklist regions (BED format)

## Usage

1. Activate the conda environment:
```bash
conda activate snakemake
```

2. Run the pipeline:
```bash
snakemake --use-conda -j <number_of_cores>
```

## Pipeline Steps

1. Quality Control
   - FastQC on raw reads

2. Read Processing
   - Adapter trimming (Trimmomatic)
   - Quality filtering

3. Alignment
   - Bowtie2 alignment
   - BAM sorting
   - Duplicate marking
   - Filtering (mapping quality, blacklist regions)

4. Peak Calling
   - MACS2 narrow peak calling

5. Differential Analysis
   - DiffBind analysis between conditions
   - Results filtered by FDR and fold change

## Output

- `results/qc/fastqc/`: FastQC HTML reports
- `results/bam/`: Processed BAM files
- `results/peaks/`: MACS2 peak files
- `results/diffbind/`: Differential binding results (CSV files)

## Notes

- The pipeline is configured to run in the "snakemake" conda environment
- Intermediate files are automatically cleaned up
- Modify thread counts in rules according to your computational resources
- Results are filtered based on FDR and fold change thresholds specified in config.yaml 