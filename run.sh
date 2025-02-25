#!/bin/bash
#SBATCH --job-name=ATAC_analysis
#SBATCH --output=logs/snakemake.out
#SBATCH --error=logs/snakemake.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Federica"
cd $WORKING_DIR || exit 1

# Create necessary directories
mkdir -p logs/cluster_logs
mkdir -p results/{qc/{fastqc,picard},trimmed,bam,peaks,diffbind}
mkdir -p tmp

# Debug information
echo "Current PATH: $PATH"
echo "Current CONDA_EXE: $CONDA_EXE"

# Unset any existing conda environment
conda deactivate 2>/dev/null || true

# Clear any existing conda initialization
unset CONDA_SHLVL
unset _CE_CONDA
unset CONDA_EXE
unset CONDA_PYTHON_EXE

# Set conda paths explicitly
export CONDA_ROOT="/opt/common/tools/ric.cosr/miniconda3"
export PATH="${CONDA_ROOT}/bin:$PATH"
export CONDA_EXE="${CONDA_ROOT}/bin/conda"

# Initialize conda
eval "$("${CONDA_ROOT}/bin/conda" shell.bash hook)"
source "${CONDA_ROOT}/etc/profile.d/conda.sh"

# Debug information after initialization
echo "Conda location: $(which conda)"
echo "Conda version: $(conda --version)"
echo "Available conda environments: "
conda env list

# Activate conda environment
conda activate snakemake

# Verify environment activation
echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python location: $(which python)"

snakemake --unlock

# Run snakemake
snakemake \
    --snakefile Snakefile \
    --executor slurm \
    --jobs 20 \
    --default-resources \
        slurm_partition=workq \
        mem_mb=8000 \
        runtime=240 \
        threads=2 \
        nodes=1 \
    --resources mem_mb=64000 \
    --set-threads bowtie2=16 mark_duplicates=4 diffbind=8 \
    --set-resources bowtie2:runtime=1440 mark_duplicates:runtime=720 diffbind:runtime=720 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix ${WORKING_DIR}/.snakemake/conda \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going