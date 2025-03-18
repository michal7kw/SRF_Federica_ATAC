#!/bin/bash
#SBATCH --job-name=create_reference
#SBATCH --output=logs/create_reference.out
#SBATCH --error=logs/create_reference.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Federica_ATAC/bowtie2_index"
cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

bowtie2-build ../hg38_reference/hg38.fa hg38
