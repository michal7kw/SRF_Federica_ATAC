library(DiffBind)
library(tidyverse)

# Get parameters from snakemake
samples <- snakemake@params$samples
comparison <- snakemake@params$comparison
fdr_threshold <- snakemake@params$fdr
lfc_threshold <- snakemake@params$lfc

# Print debugging information about input parameters
print("Input parameters:")
print("Samples:")
str(samples)
print("Comparison:")
str(comparison)

# Get relevant samples for this comparison
treatment <- comparison$treatment
control <- comparison$control
relevant_conditions <- c(treatment, control)

# Get all sample names
all_samples <- names(samples)
relevant_samples <- all_samples[sapply(all_samples, function(x) samples[[x]]$condition %in% relevant_conditions)]

print("Relevant samples for comparison:")
print(relevant_samples)

# Create sample metadata
sample_metadata <- data.frame(
    SampleID = relevant_samples,
    Condition = sapply(relevant_samples, function(x) samples[[x]]$condition),
    Replicate = seq_along(relevant_samples),
    bamReads = file.path(getwd(), sprintf("results/bam/%s.filtered.bam", relevant_samples)),
    Peaks = file.path(getwd(), sprintf("results/peaks/%s_peaks.narrowPeak", relevant_samples)),
    PeakCaller = "narrow",
    stringsAsFactors = FALSE
)

# Print debugging information
print("Sample metadata:")
print(sample_metadata)

# Check if peak files exist and have content
print("Checking peak files:")
for(peak_file in sample_metadata$Peaks) {
    if(!file.exists(peak_file)) {
        stop(sprintf("Peak file does not exist: %s", peak_file))
    }
    peak_content <- try(read.table(peak_file, header=FALSE), silent=TRUE)
    if(inherits(peak_content, "try-error")) {
        stop(sprintf("Error reading peak file %s: %s", peak_file, attr(peak_content, "condition")$message))
    }
    if(nrow(peak_content) == 0) {
        stop(sprintf("Peak file %s is empty", peak_file))
    }
    print(sprintf("%s: %d peaks", basename(peak_file), nrow(peak_content)))
}

# Check BAM files
print("Checking BAM files:")
for(bam_file in sample_metadata$bamReads) {
    if(!file.exists(bam_file)) {
        stop(sprintf("BAM file does not exist: %s", bam_file))
    }
    print(sprintf("Found BAM file: %s", basename(bam_file)))
}

# Create DiffBind object
print("Creating DiffBind object...")
dba <- dba(sampleSheet=sample_metadata)
print("DiffBind object created successfully")
print(dba)

# Count reads
print("Counting reads...")
dba <- dba.count(dba)
print("Read counting complete")
print(dba)

# Set up contrasts
print("Setting up contrasts...")
# Ensure group assignments are correct
group1_samples <- which(dba$samples$Condition == comparison$treatment)
group2_samples <- which(dba$samples$Condition == comparison$control)

print("Group assignments:")
print(paste("Treatment group (", comparison$treatment, "):", paste(group1_samples, collapse=", ")))
print(paste("Control group (", comparison$control, "):", paste(group2_samples, collapse=", ")))

if(length(group1_samples) == 0 || length(group2_samples) == 0) {
    stop("Error: One or both groups have no samples")
}

# Create contrast
dba <- dba.contrast(dba,
                    group1=group1_samples,
                    group2=group2_samples,
                    name1=comparison$treatment,
                    name2=comparison$control)

print("Contrast created:")
print(dba$contrasts[[1]])

# Perform differential analysis
print("Running differential analysis...")
dba <- dba.analyze(dba, method=DBA_DESEQ2)

# Extract results
print("Extracting results...")
res <- dba.report(dba)
print(paste("Found", length(res), "differential peaks"))

# Convert to data frame and write output
print("Writing results...")
res_df <- as.data.frame(res)
write.csv(res_df, snakemake@output[[1]], row.names=FALSE)

print("Analysis complete") 