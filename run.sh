#!/bin/bash
#SBATCH --job-name=F_ATAC_analysis
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
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Federica_ATAC"
cd $WORKING_DIR || exit 1

# Create all necessary directories
mkdir -p logs/cluster_logs
mkdir -p logs/filter_bam
mkdir -p logs/create_bigwig
mkdir -p logs/tss_enrichment
mkdir -p logs/annotate_peaks
mkdir -p results/{qc/{fastqc,picard,fragment_sizes,tss_enrichment,multiqc},trimmed,bam,peaks,diffbind,tracks,annotation,visualization}
mkdir -p tmp
mkdir -p envs

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

# Unlock the working directory if needed
snakemake --unlock

# Verify required input files exist
echo "Checking for required input files..."
for SAMPLE in $(grep -v '^#' config.yaml | grep -o '"[^"]*": {' | sed 's/": {//' | sed 's/"//g'); do
    BAM_FILE="results/bam/${SAMPLE}.markdup.bam"
    if [ ! -f "$BAM_FILE" ]; then
        echo "Warning: $BAM_FILE does not exist. Pipeline may fail."
    fi
done

if [ ! -f "$(grep 'blacklist' config.yaml | awk '{print $2}' | tr -d '"')" ]; then
    echo "Warning: Blacklist file not found. Pipeline may fail."
fi

if [ ! -f "$(grep 'tss_bed' config.yaml | awk '{print $2}' | tr -d '"')" ]; then
    echo "Warning: TSS bed file not found. Pipeline may fail."
fi

# Create conda environments if they don't exist
echo "===== Creating conda environments ====="
snakemake --use-conda --conda-create-envs-only --conda-frontend conda \
    --conda-prefix ${WORKING_DIR}/.snakemake/conda

# Run snakemake in multiple phases to ensure proper order
echo "===== PHASE 1: Filtering BAM files ====="
snakemake \
    --snakefile Snakefile \
    --executor slurm \
    --jobs 12 \
    --default-resources \
        slurm_partition=workq \
        mem_mb=8000 \
        runtime=240 \
        threads=2 \
        nodes=1 \
    --resources mem_mb=64000 \
    --set-threads bowtie2=16 mark_duplicates=4 filter_bam=8 \
    --set-resources bowtie2:runtime=1440 mark_duplicates:runtime=720 filter_bam:runtime=240 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix ${WORKING_DIR}/.snakemake/conda \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --until filter_bam

echo "===== PHASE 2: Creating bigWig files ====="
snakemake \
    --snakefile Snakefile \
    --executor slurm \
    --jobs 12 \
    --default-resources \
        slurm_partition=workq \
        mem_mb=8000 \
        runtime=240 \
        threads=2 \
        nodes=1 \
    --resources mem_mb=64000 \
    --set-threads create_bigwig=8 \
    --set-resources create_bigwig:runtime=240 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix ${WORKING_DIR}/.snakemake/conda \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --until create_bigwig

echo "===== PHASE 3: Running analyses that depend on bigWig files ====="
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
    --set-threads diffbind=8 \
    --set-resources diffbind:runtime=720 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix ${WORKING_DIR}/.snakemake/conda \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going

echo "===== Pipeline execution complete ====="