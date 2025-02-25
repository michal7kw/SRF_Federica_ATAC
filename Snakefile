import os
from snakemake.utils import min_version

# Set minimum snakemake version
min_version("6.0")

# Load configuration
configfile: "config.yaml"

# Helper functions
def get_final_output():
    final_output = []
    # QC reports
    final_output.extend(expand("results/qc/fastqc/{sample}_{read}_fastqc.html", 
                             sample=config["samples"], 
                             read=["r1", "r2"]))
    # BAM files
    final_output.extend(expand("results/bam/{sample}.filtered.bam", 
                             sample=config["samples"]))
    # MACS2 peaks
    final_output.extend(expand("results/peaks/{sample}_peaks.narrowPeak", 
                             sample=config["samples"]))
    # DiffBind results
    final_output.extend(expand("results/diffbind/{comparison}/differential_peaks.csv", 
                             comparison=[comp["name"] for comp in config["comparisons"]]))
    return final_output

# Rules
rule all:
    input:
        get_final_output()

# FastQC - Light on resources, can run in parallel
rule fastqc:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        html_r1 = "results/qc/fastqc/{sample}_r1_fastqc.html",
        html_r2 = "results/qc/fastqc/{sample}_r2_fastqc.html",
        zip_r1 = "results/qc/fastqc/{sample}_r1_fastqc.zip",
        zip_r2 = "results/qc/fastqc/{sample}_r2_fastqc.zip"
    resources:
        mem_mb=4000,
        runtime=60
    threads: 2
    conda:
        "snakemake"
    shell:
        """
        # Run FastQC
        mkdir -p results/qc/fastqc
        fastqc {input.r1} {input.r2} -o results/qc/fastqc/ -t {threads}
        
        # Rename the output files to match Snakemake's expected names
        for f in results/qc/fastqc/*_R1_001_fastqc.*; do
            if [ -f "$f" ]; then
                newname=$(echo $f | sed 's/_R1_001_fastqc/_r1_fastqc/')
                mv "$f" "$newname"
            fi
        done
        
        for f in results/qc/fastqc/*_R2_001_fastqc.*; do
            if [ -f "$f" ]; then
                newname=$(echo $f | sed 's/_R2_001_fastqc/_r2_fastqc/')
                mv "$f" "$newname"
            fi
        done
        """

# Trimmomatic - Memory intensive for large files
rule trimmomatic:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        r1_paired = "results/trimmed/{sample}_r1_paired.fastq.gz",
        r2_paired = "results/trimmed/{sample}_r2_paired.fastq.gz",
        r1_unpaired = "results/trimmed/{sample}_r1_unpaired.fastq.gz",
        r2_unpaired = "results/trimmed/{sample}_r2_unpaired.fastq.gz"
    params:
        adapters = config["params"]["trimmomatic"]["adapters"],
        leading = config["params"]["trimmomatic"]["leading"],
        trailing = config["params"]["trimmomatic"]["trailing"],
        slidingwindow = config["params"]["trimmomatic"]["slidingwindow"],
        minlen = config["params"]["trimmomatic"]["minlen"]
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    conda:
        "snakemake"
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 \
            {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            SLIDINGWINDOW:{params.slidingwindow} \
            MINLEN:{params.minlen}
        """

# Bowtie2 alignment - CPU and memory intensive
rule bowtie2:
    input:
        r1 = "results/trimmed/{sample}_r1_paired.fastq.gz",
        r2 = "results/trimmed/{sample}_r2_paired.fastq.gz"
    output:
        bam = temp("results/bam/{sample}.raw.bam")
    params:
        index = config["genome"]["bowtie2_index"]
    threads: 16
    resources:
        mem_mb=32000,
        runtime=1440
    conda:
        "snakemake"
    shell:
        """
        # Create temporary directory
        mkdir -p tmp/{wildcards.sample}
        
        # Create read group string
        RG="@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA"
        
        # Run bowtie2 with read group information and proper BAM creation
        bowtie2 -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} --rg-id {wildcards.sample} \
            --rg "SM:{wildcards.sample}" --rg "LB:{wildcards.sample}" --rg "PL:ILLUMINA" 2> tmp/{wildcards.sample}/align.log | \
            samtools view -bS - | \
            samtools view -h -b > {output.bam}
            
        # Verify the BAM file is valid
        samtools quickcheck {output.bam}
        
        # Clean up
        rm -rf tmp/{wildcards.sample}
        """

# Sort BAM - I/O intensive, benefits from multiple threads
rule sort_bam:
    input:
        "results/bam/{sample}.raw.bam"
    output:
        bam = temp("results/bam/{sample}.sorted.bam"),
        bai = temp("results/bam/{sample}.sorted.bam.bai")
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    conda:
        "snakemake"
    shell:
        """
        # Create temporary directory
        mkdir -p tmp/{wildcards.sample}
        
        # Verify input BAM is valid
        samtools quickcheck {input}
        
        # Sort BAM file
        samtools sort -@ {threads} \
            -m 2G \
            -T tmp/{wildcards.sample}/sort \
            -o {output.bam} {input}
            
        # Index BAM file
        samtools index {output.bam}
        
        # Verify read groups are present
        samtools view -H {output.bam} | grep "^@RG"
        
        # Clean up
        rm -rf tmp/{wildcards.sample}
        """

# Mark duplicates - Memory intensive
rule mark_duplicates:
    input:
        "results/bam/{sample}.sorted.bam"
    output:
        bam = "results/bam/{sample}.markdup.bam",
        metrics = "results/qc/picard/{sample}.markdup_metrics.txt"
    resources:
        mem_mb=32000,
        runtime=720
    threads: 4
    conda:
        "snakemake"
    shell:
        """
        # Create temporary directory
        mkdir -p tmp/{wildcards.sample}

        # Add read groups if missing
        if ! samtools view -H {input} | grep -q "^@RG"; then
            # Create read group string
            RG="@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA"
            
            # Add read group to header
            samtools addreplacerg -r "$RG" -o tmp/{wildcards.sample}/with_rg.bam {input}
            INPUT_BAM="tmp/{wildcards.sample}/with_rg.bam"
        else
            INPUT_BAM="{input}"
        fi

        # Run Picard MarkDuplicates
        picard -Xmx{resources.mem_mb}m MarkDuplicates \
            INPUT=$INPUT_BAM \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            REMOVE_DUPLICATES=false \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR=tmp/{wildcards.sample}

        # Clean up
        rm -rf tmp/{wildcards.sample}
        """

# Filter BAM - CPU intensive for large files
rule filter_bam:
    input:
        bam = "results/bam/{sample}.markdup.bam",
        blacklist = config["genome"]["blacklist"]
    output:
        "results/bam/{sample}.filtered.bam"
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    conda:
        "snakemake"
    shell:
        """
        samtools view -@ {threads} -b -F 1804 -q 30 {input.bam} | \
        bedtools intersect -v -a stdin -b {input.blacklist} > {output}
        samtools index -@ {threads} {output}
        """

# Call peaks with MACS2 - Memory intensive
rule macs2_callpeak:
    input:
        "results/bam/{sample}.filtered.bam"
    output:
        peaks = "results/peaks/{sample}_peaks.narrowPeak",
        xls = "results/peaks/{sample}_peaks.xls",
        summits = "results/peaks/{sample}_summits.bed"
    params:
        format = config["params"]["macs2"]["format"],
        genome_size = "hs",
        qvalue = config["params"]["macs2"]["qvalue"],
        name = "{sample}"
    resources:
        mem_mb=32000,  # Increased memory
        runtime=240
    threads: 4
    conda:
        "snakemake"
    shell:
        """
        # Create output and temporary directories
        mkdir -p results/peaks
        mkdir -p tmp/{wildcards.sample}
        
        # Run MACS2 with increased buffer size and explicit temp directory
        macs2 callpeak \
            -t {input} \
            -f {params.format} \
            -g {params.genome_size} \
            -n {params.name} \
            --outdir results/peaks \
            -q {params.qvalue} \
            --tempdir tmp/{wildcards.sample} \
            --buffer-size 1000000 || true  # Continue even if there are non-critical errors
            
        # Verify output files exist and are not empty
        if [ ! -s {output.peaks} ]; then
            echo "Error: Peak file {output.peaks} is empty or missing"
            exit 1
        fi
        if [ ! -s {output.xls} ]; then
            echo "Error: XLS file {output.xls} is empty or missing"
            exit 1
        fi
        if [ ! -s {output.summits} ]; then
            echo "Error: Summits file {output.summits} is empty or missing"
            exit 1
        fi
        """

# DiffBind analysis - Very memory intensive
rule diffbind:
    input:
        peaks = lambda wildcards: expand(
            "results/peaks/{sample}_peaks.narrowPeak",
            sample=[s for s in config["samples"] if config["samples"][s]["condition"] in 
                   [next(c["treatment"] for c in config["comparisons"] if c["name"] == wildcards.comparison),
                    next(c["control"] for c in config["comparisons"] if c["name"] == wildcards.comparison)]]
        ),
        bams = lambda wildcards: expand(
            "results/bam/{sample}.filtered.bam",
            sample=[s for s in config["samples"] if config["samples"][s]["condition"] in 
                   [next(c["treatment"] for c in config["comparisons"] if c["name"] == wildcards.comparison),
                    next(c["control"] for c in config["comparisons"] if c["name"] == wildcards.comparison)]]
        )
    output:
        "results/diffbind/{comparison}/differential_peaks.csv"
    params:
        samples = config["samples"],
        comparison = lambda wildcards: next(comp for comp in config["comparisons"] if comp["name"] == wildcards.comparison),
        fdr = config["params"]["diffbind"]["fdr"],
        lfc = config["params"]["diffbind"]["lfc"]
    resources:
        mem_mb=64000,
        runtime=480
    threads: 8
    conda:
        "snakemake"
    script:
        "scripts/run_diffbind.R" 