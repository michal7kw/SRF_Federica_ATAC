#!/bin/bash
#SBATCH --job-name=get_tss_from_genecode
#SBATCH --output=logs/get_tss_from_genecode.out
#SBATCH --error=logs/get_tss_from_genecode.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Federica_ATAC/resources"
cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

GetTss -d gencode -g ../data/gencode.v47.basic.annotation.gtf -t gencode_tss.bed


# https://github.com/junjunlab/GetTss
