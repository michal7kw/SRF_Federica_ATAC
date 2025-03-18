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
    # Add MultiQC report to final output
    final_output.append("results/qc/multiqc/multiqc_report.html")
    # Add bigWig tracks to final output
    final_output.extend(expand("results/tracks/{sample}.bw", 
                             sample=config["samples"]))
    # Add fragment size QC
    final_output.extend(expand("results/qc/fragment_sizes/{sample}.fragsize.pdf", 
                             sample=config["samples"]))
    # Add TSS enrichment QC
    final_output.extend(expand("results/qc/tss_enrichment/{sample}.tss_enrichment.pdf", 
                             sample=config["samples"]))
    # Add peak annotation output
    final_output.extend(expand("results/annotation/{sample}_annotated_peaks.tsv", 
                             sample=config["samples"]))
    # Add peak visualization output
    final_output.extend(expand("results/visualization/{sample}_peak_heatmap.pdf", 
                             sample=config["samples"]))
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
        # Create temporary directory and output directory
        mkdir -p tmp/{wildcards.sample}
        mkdir -p results/bam
        
        # Verify input files exist and are not empty
        echo "Checking input files..."
        ls -lh {input.r1} {input.r2}
        
        if [ ! -s {input.r1} ]; then echo "Error: {input.r1} does not exist or is empty" && exit 1; fi
        if [ ! -s {input.r2} ]; then echo "Error: {input.r2} does not exist or is empty" && exit 1; fi
        
        # Check if index exists and list its contents
        echo "Checking bowtie2 index..."
        if [ ! -d $(dirname {params.index}) ]; then 
            echo "Error: Bowtie2 index directory not found at $(dirname {params.index})" && exit 1
        else
            echo "Index directory contents:"
            ls -la $(dirname {params.index})
        fi
        
        # Check for index files specifically
        echo "Looking for .bt2 files in index directory..."
        find $(dirname {params.index}) -name "*.bt2" -o -name "*.bt2l"
        
        echo "Starting bowtie2 alignment for {wildcards.sample}..."
        
        # Step 1: Run bowtie2 to create SAM file with better error handling
        echo "Running bowtie2 command..."
        
        # First check bowtie2 version
        echo "Bowtie2 version:"
        bowtie2 --version
        
        # Run with more verbose output
        bowtie2 -p {threads} \
            -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -X 2000 \
            --very-sensitive \
            --rg-id "{wildcards.sample}" \
            --rg "SM:{wildcards.sample}" \
            --rg "LB:{wildcards.sample}" \
            --rg "PL:ILLUMINA" \
            -S tmp/{wildcards.sample}/aligned.sam \
            --no-unal \
            --time \
            2> tmp/{wildcards.sample}/align.log || {{ 
                echo "Bowtie2 alignment failed. Error log:" 
                cat tmp/{wildcards.sample}/align.log
                exit 1
            }}
        
        # Check if SAM file was created successfully
        if [ ! -s tmp/{wildcards.sample}/aligned.sam ]; then
            echo "Error: SAM file not created or empty"
            echo "Bowtie2 log contents:"
            cat tmp/{wildcards.sample}/align.log
            exit 1
        fi
        
        # Step 2: Convert SAM to BAM
        echo "Converting SAM to BAM..."
        samtools view -@ 4 -bS tmp/{wildcards.sample}/aligned.sam > tmp/{wildcards.sample}/unsorted.bam || {{
            echo "Error: SAM to BAM conversion failed"
            exit 1
        }}
        
        # Check if BAM conversion was successful
        if [ ! -s tmp/{wildcards.sample}/unsorted.bam ]; then
            echo "Error: BAM conversion failed"
            exit 1
        fi
        
        # Step 3: Sort BAM file
        echo "Sorting BAM file..."
        samtools sort -@ 4 -m 4G -T tmp/{wildcards.sample}/sort \
            -o {output.bam} tmp/{wildcards.sample}/unsorted.bam || {{
            echo "Error: BAM sorting failed"
            exit 1
        }}
        
        # Verify the BAM file was created successfully
        if [ ! -s {output.bam} ]; then
            echo "Error: Sorted BAM file not created or empty"
            exit 1
        fi
        
        # Check the BAM file
        echo "Validating BAM file..."
        samtools quickcheck {output.bam}
        if [ $? -ne 0 ]; then
            echo "Error: BAM file failed validation"
            samtools view -H {output.bam}
            exit 1
        fi
        
        # Index the BAM file to make sure it's valid
        echo "Indexing BAM file..."
        samtools index {output.bam}
        
        # Clean up temporary files
        rm -f tmp/{wildcards.sample}/aligned.sam
        rm -f tmp/{wildcards.sample}/unsorted.bam
        rm -rf tmp/{wildcards.sample}
        
        echo "Bowtie2 alignment for {wildcards.sample} completed successfully"
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
        bam = "results/bam/{sample}.filtered.bam",
        bai = "results/bam/{sample}.filtered.bam.bai"
    params:
        # Mitochondrial contig name may vary by genome assembly
        mito_chrom = config["genome"]["mito_chrom"],  # Usually "chrM" or "MT"
        # Use sample-specific temporary files
        tmp_filtered = "tmp/{sample}.filtered.tmp.bam",
        tmp_nomito = "tmp/{sample}.nomito.tmp.bam"
    log:
        "logs/filter_bam/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    conda:
        "snakemake"
    shell:
        """
        # Create necessary directories
        mkdir -p tmp
        mkdir -p logs/filter_bam
        
        # Step 1: Filter BAM for quality metrics
        echo "Filtering for quality metrics..." > {log}
        samtools view -@ {threads} -b -F 1804 -f 2 -q 30 \
            -o {params.tmp_filtered} {input.bam} 2>> {log}
        
        # Step 2: Remove mitochondrial reads (fixed approach)
        echo "Removing mitochondrial reads..." >> {log}
        samtools view -h {params.tmp_filtered} 2>> {log} | \
        grep -v "^@" | grep -v "{params.mito_chrom}\t" | \
        cat <(samtools view -H {params.tmp_filtered} 2>> {log}) - | \
        samtools view -@ {threads} -b -o {params.tmp_nomito} 2>> {log}
        
        # Step 3: Remove blacklisted regions
        echo "Removing blacklisted regions..." >> {log}
        bedtools intersect -v -a {params.tmp_nomito} -b {input.blacklist} > {output.bam} 2>> {log}
        
        # Step 4: Index the filtered BAM
        echo "Indexing filtered BAM..." >> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        
        # Clean up temporary files
        echo "Cleaning up..." >> {log}
        rm -f {params.tmp_filtered} {params.tmp_nomito}
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
        format = "BAMPE",  # Paired-end format for ATAC-seq
        genome_size = "hs",
        qvalue = config["params"]["macs2"]["qvalue"],
        name = "{sample}",
        # ATAC-seq specific parameters
        shift = "-75",
        extsize = "150",
        nomodel = True
    resources:
        mem_mb=32000,
        runtime=240
    threads: 4
    conda:
        "snakemake"
    shell:
        """
        # Create output and temporary directories
        mkdir -p results/peaks
        mkdir -p tmp/{wildcards.sample}
        
        # Run MACS2 with ATAC-seq specific parameters
        macs2 callpeak \
            -t {input} \
            -f {params.format} \
            -g {params.genome_size} \
            -n {params.name} \
            --outdir results/peaks \
            -q {params.qvalue} \
            --shift {params.shift} \
            --extsize {params.extsize} \
            --nomodel \
            --keep-dup all \
            --tempdir tmp/{wildcards.sample} \
            --buffer-size 1000000 || true
            
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

# MultiQC - Aggregate QC reports
rule multiqc:
    input:
        fastqc = expand("results/qc/fastqc/{sample}_{read}_fastqc.zip", 
                      sample=config["samples"], read=["r1", "r2"]),
        markdup = expand("results/qc/picard/{sample}.markdup_metrics.txt", 
                       sample=config["samples"]),
        fragment_sizes = expand("results/qc/fragment_sizes/{sample}.fragsize.txt", 
                              sample=config["samples"]),
        tss_enrichment = expand("results/qc/tss_enrichment/{sample}.tss_matrix.gz", 
                              sample=config["samples"])
    output:
        report = "results/qc/multiqc/multiqc_report.html"
    params:
        outdir = "results/qc/multiqc"
    resources:
        mem_mb=8000,
        runtime=60
    threads: 1
    conda:
        "snakemake"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run MultiQC
        multiqc \
            results/qc/fastqc \
            results/qc/picard \
            results/qc/fragment_sizes \
            results/qc/tss_enrichment \
            -o {params.outdir} \
            -f \
            -v
        """

# Generate bigWig files for genome browser visualization
rule create_bigwig:
    input:
        bam = "results/bam/{sample}.filtered.bam",
        bai = "results/bam/{sample}.filtered.bam.bai"
    output:
        bigwig = "results/tracks/{sample}.bw"
    params:
        genome_size = config["genome"]["chrom_sizes"],
        normalize = "RPKM"  # Options: CPM, RPKM, BPM, RPGC, None
    log:
        "logs/create_bigwig/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    conda:
        "deeptools_env"
    shell:
        """
        # Create output directory
        mkdir -p results/tracks
        mkdir -p logs/create_bigwig
        
        # Verify input BAM file
        echo "Checking input BAM file..." > {log}
        samtools quickcheck {input.bam} || {{ echo "ERROR: BAM file {input.bam} failed validation" >> {log}; exit 1; }}
        
        # Generate normalized bigWig file with ATAC-seq specific parameters
        echo "Running bamCoverage..." >> {log}
        bamCoverage -b {input.bam} -o {output.bigwig} \
            --binSize 10 \
            --normalizeUsing {params.normalize} \
            --effectiveGenomeSize $(cat {params.genome_size} | awk '{{sum+=$2}} END {{print sum}}') \
            --numberOfProcessors {threads} \
            --ignoreDuplicates \
            --minMappingQuality 30 \
            --centerReads \
            2>> {log}
            
        # Verify output bigWig file
        echo "Checking output bigWig file..." >> {log}
        if [ ! -s {output.bigwig} ]; then
            echo "ERROR: BigWig file {output.bigwig} was not created or is empty" >> {log}
            exit 1
        fi
        
        echo "BigWig file generated successfully" >> {log}
        """

# Fragment size analysis - important QC metric for ATAC-seq
rule fragment_size_analysis:
    input:
        bam = "results/bam/{sample}.filtered.bam"
    output:
        txt = "results/qc/fragment_sizes/{sample}.fragsize.txt",
        pdf = "results/qc/fragment_sizes/{sample}.fragsize.pdf"
    resources:
        mem_mb=8000,
        runtime=120
    threads: 4
    conda:
        "snakemake"
    shell:
        """
        # Create output directory
        mkdir -p results/qc/fragment_sizes
        
        # Extract fragment lengths and generate histogram
        samtools view -f 2 {input.bam} | \
        awk '{{if ($9>0) print $9}}' | \
        sort -n | uniq -c > {output.txt}
        
        # Plot fragment size distribution using R
        Rscript -e '
        data <- read.table("{output.txt}", header=FALSE)
        colnames(data) <- c("Count", "Size")
        pdf("{output.pdf}", width=8, height=6)
        plot(data$Size, data$Count, type="l", col="blue", 
             main="Fragment Size Distribution for {wildcards.sample}",
             xlab="Fragment Size (bp)", ylab="Count")
        abline(v=c(50, 100, 150, 200, 250), lty=2, col="gray")
        dev.off()
        '
        """

# TSS enrichment analysis - key QC metric for ATAC-seq
rule tss_enrichment:
    input:
        bam = "results/bam/{sample}.filtered.bam",
        bai = "results/bam/{sample}.filtered.bam.bai",
        bigwig = "results/tracks/{sample}.bw"  # Add explicit dependency on bigwig
    output:
        matrix = "results/qc/tss_enrichment/{sample}.tss_matrix.gz",
        plot = "results/qc/tss_enrichment/{sample}.tss_enrichment.pdf"
    params:
        tss_bed = config["genome"]["tss_bed"],
        before = 2000,  # bp upstream of TSS
        after = 2000    # bp downstream of TSS
    log:
        "logs/tss_enrichment/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    conda:
        "deeptools_env"
    shell:
        """
        # Create output directory
        mkdir -p results/qc/tss_enrichment
        mkdir -p logs/tss_enrichment
        
        # Check if bigWig file exists and is not empty
        if [ ! -s {input.bigwig} ]; then
            echo "ERROR: BigWig file {input.bigwig} does not exist or is empty" > {log}
            exit 1
        fi
        
        # Calculate signal matrix around TSS
        echo "Running computeMatrix..." >> {log}
        computeMatrix reference-point \
            --referencePoint TSS \
            -b {params.before} -a {params.after} \
            -R {params.tss_bed} \
            -S {input.bigwig} \
            --skipZeros \
            --numberOfProcessors {threads} \
            -o {output.matrix} \
            2>> {log}
            
        # Plot the TSS enrichment profile
        echo "Running plotProfile..." >> {log}
        plotProfile \
            -m {output.matrix} \
            --perGroup \
            --refPointLabel "TSS" \
            -out {output.plot} \
            --plotTitle "TSS Enrichment - {wildcards.sample}" \
            --yAxisLabel "Signal" \
            --regionsLabel "TSS" \
            2>> {log}
        
        echo "TSS enrichment analysis completed successfully" >> {log}
        """

# Peak annotation with ChIPseeker
rule annotate_peaks:
    input:
        peaks = "results/peaks/{sample}_peaks.narrowPeak"
    output:
        annotated = "results/annotation/{sample}_annotated_peaks.tsv",
        plot = "results/annotation/{sample}_annotation_plot.pdf"
    params:
        txdb = config["genome"]["txdb"],
        genome = config["genome"]["name"]
    log:
        "logs/annotate_peaks/{sample}.log"
    resources:
        mem_mb=16000,
        runtime=120
    threads: 2
    conda:
        "envs/chipseeker.yaml"  # Create this environment file
    script:
        "scripts/annotate_peaks.R"

# Generate peak heatmap
rule peak_heatmap:
    input:
        bw = "results/tracks/{sample}.bw",
        peaks = "results/peaks/{sample}_peaks.narrowPeak"
    output:
        matrix = "results/visualization/{sample}_peak_matrix.gz",
        heatmap = "results/visualization/{sample}_peak_heatmap.pdf",
        profile = "results/visualization/{sample}_peak_profile.pdf"
    params:
        window_size = 1000  # Changed from 'extend' to 'window_size'
    resources:
        mem_mb=16000,
        runtime=240
    threads: 8
    conda:
        "deeptools_env"
    shell:
        """
        # Create output directory
        mkdir -p results/visualization
        
        # Convert narrowPeak to BED format for deepTools
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4}}' {input.peaks} > tmp_{wildcards.sample}.bed
        
        # Compute matrix around peak centers
        computeMatrix reference-point \
            --referencePoint center \
            -b {params.window_size} -a {params.window_size} \
            -R tmp_{wildcards.sample}.bed \
            -S {input.bw} \
            --skipZeros \
            --numberOfProcessors {threads} \
            -o {output.matrix}
            
        # Generate heatmap
        plotHeatmap \
            -m {output.matrix} \
            -out {output.heatmap} \
            --colorMap YlOrRd \
            --whatToShow 'heatmap and colorbar' \
            --heatmapHeight 15 \
            --heatmapWidth 4 \
            --plotTitle "{wildcards.sample} Peak Enrichment"
            
        # Generate profile plot
        plotProfile \
            -m {output.matrix} \
            -out {output.profile} \
            --plotTitle "{wildcards.sample} Average Peak Profile" \
            --perGroup
            
        # Clean up
        rm -f tmp_{wildcards.sample}.bed
        """ 