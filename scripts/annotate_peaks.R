# Load required libraries with error handling
suppressPackageStartupMessages({
  if (!require("ChIPseeker")) {
    stop("ChIPseeker package is required but not installed")
  }
  if (!require("GenomicFeatures")) {
    stop("GenomicFeatures package is required but not installed")
  }
})

# Get parameters from Snakemake
peaks_file <- snakemake@input[["peaks"]]
output_tsv <- snakemake@output[["annotated"]]
output_plot <- snakemake@output[["plot"]]
genome <- snakemake@params[["genome"]]
txdb_name <- snakemake@params[["txdb"]]

# Load appropriate annotation packages based on genome
tryCatch({
  if (txdb_name == "hg38") {
    if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene")) {
      stop("TxDb.Hsapiens.UCSC.hg38.knownGene package is required but not installed")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    if (!require("org.Hs.eg.db")) {
      stop("org.Hs.eg.db package is required but not installed")
    }
    orgdb <- org.Hs.eg.db
  } else if (txdb_name == "hg19") {
    if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene")) {
      stop("TxDb.Hsapiens.UCSC.hg19.knownGene package is required but not installed")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    
    if (!require("org.Hs.eg.db")) {
      stop("org.Hs.eg.db package is required but not installed")
    }
    orgdb <- org.Hs.eg.db
  } else if (txdb_name == "mm10") {
    if (!require("TxDb.Mmusculus.UCSC.mm10.knownGene")) {
      stop("TxDb.Mmusculus.UCSC.mm10.knownGene package is required but not installed")
    }
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    
    if (!require("org.Mm.eg.db")) {
      stop("org.Mm.eg.db package is required but not installed")
    }
    orgdb <- org.Mm.eg.db
  } else {
    stop(paste0("Unsupported genome: ", txdb_name))
  }
}, error = function(e) {
  cat("Error loading annotation packages:", conditionMessage(e), "\n")
  quit(status = 1)
})

# Read and validate peaks file
tryCatch({
  cat(">> Reading peaks file...\t", format(Sys.time(), "%Y-%m-%d %I:%M:%S %p"), "\n")
  # Check if file exists and is not empty
  if (!file.exists(peaks_file) || file.size(peaks_file) == 0) {
    stop(paste0("Peak file does not exist or is empty: ", peaks_file))
  }
  
  # Validate file format
  con <- file(peaks_file, "r")
  first_line <- readLines(con, n = 1)
  close(con)
  
  # Check if it's a valid narrowPeak file (should have at least 6 columns)
  if (length(strsplit(first_line, "\t")[[1]]) < 6) {
    stop("Invalid narrowPeak format. File should be tab-delimited with at least 6 columns.")
  }
  
  # Read peaks with error handling
  peaks <- readPeakFile(peaks_file)
  
  # Verify peaks were loaded
  if (length(peaks) == 0) {
    stop("No peaks were loaded from the input file.")
  }
  
  cat(">> preparing features information...\t", format(Sys.time(), "%Y-%m-%d %I:%M:%S %p"), "\n")
  
  # Add chromosome lengths to prevent out of bounds errors
  seqlevels_peaks <- seqlevels(peaks)
  seqinfo_txdb <- seqinfo(txdb)
  
  # Make sure we have compatible chromosome naming
  if (!any(seqlevels_peaks %in% seqlevels(seqinfo_txdb))) {
    # Try to standardize chromosome names (e.g., convert "chr1" to "1" or vice versa)
    if (all(grepl("^chr", seqlevels_peaks))) {
      # If peaks use "chr1" but txdb uses "1"
      seqlevels(peaks) <- sub("^chr", "", seqlevels_peaks)
    } else {
      # If peaks use "1" but txdb uses "chr1"
      seqlevels(peaks) <- paste0("chr", seqlevels_peaks)
    }
    
    # Check again after conversion
    if (!any(seqlevels(peaks) %in% seqlevels(seqinfo_txdb))) {
      stop("Chromosome naming in peak file is incompatible with annotation database")
    }
  }
  
  cat(">> identifying nearest features...\t", format(Sys.time(), "%Y-%m-%d %I:%M:%S %p"), "\n")
  
  # Annotate peaks
  peakAnno <- annotatePeak(peaks, 
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         annoDb = orgdb)
  
  cat(">> writing annotation results...\t", format(Sys.time(), "%Y-%m-%d %I:%M:%S %p"), "\n")
  
  # Export annotation to TSV file
  anno_df <- as.data.frame(peakAnno)
  write.table(anno_df, file = output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Generate annotation plots
  cat(">> generating annotation plots...\t", format(Sys.time(), "%Y-%m-%d %I:%M:%S %p"), "\n")
  pdf(output_plot, width = 12, height = 10)
  
  # Plot annotation statistics
  plotAnnoPie(peakAnno)
  
  # Plot distribution of peaks relative to TSS
  plotDistToTSS(peakAnno, title = "Distribution of peaks relative to TSS")
  
  dev.off()
  
  cat(">> annotation complete!\t", format(Sys.time(), "%Y-%m-%d %I:%M:%S %p"), "\n")
  
}, error = function(e) {
  cat("Error in peak annotation:", conditionMessage(e), "\n")
  quit(status = 1)
}) 