# ============================================================
# 09_align_reads.R
#
# Align trimmed paired-end RNA-seq reads to GRCh38 using
# Rsubread::subjunc().
# ============================================================

setwd("C:/Users/benma/RNA-seq")

library(Rsubread)
library(tidyverse)

source("config/config.R")


# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

fastq_dir <- file.path(
  "data",
  "raw",
  "fastq",
  "trimmed"
)

reference_dir <- file.path(
  "reference",
  "genome"
)

bam_dir <- file.path(
  "data",
  "processed",
  "bam"
)

log_dir <- file.path(
  "results",
  "logs"
)


dir.create(
  bam_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Absolute paths
# ------------------------------------------------------------

fastq_dir_abs <- normalizePath(
  fastq_dir,
  winslash = "/",
  mustWork = TRUE
)

bam_dir_abs <- normalizePath(
  bam_dir,
  winslash = "/",
  mustWork = TRUE
)

index_base <- normalizePath(
  file.path(
    reference_dir,
    "GRCh38_rsubread"
  ),
  winslash = "/",
  mustWork = FALSE
)


# ------------------------------------------------------------
# Samples
# ------------------------------------------------------------

samples <- TEST_SAMPLES

cat("Samples to align:\n")
print(samples)


# ------------------------------------------------------------
# Check index
# ------------------------------------------------------------

index_files <- list.files(
  dirname(index_base),
  pattern = "^GRCh38_rsubread",
  full.names = TRUE
)

if (length(index_files) == 0) {
  
  stop(
    "Rsubread index was not found."
  )
  
}

cat("\nRsubread index found:\n")
print(index_files)


# ------------------------------------------------------------
# Align samples
# ------------------------------------------------------------

for (sample in samples) {
  
  cat("\n========================================\n")
  cat("Aligning sample:", sample, "\n")
  cat("========================================\n")
  
  
  # ----------------------------------------------------------
  # Input FASTQ
  # ----------------------------------------------------------
  
  read1 <- file.path(
    fastq_dir_abs,
    paste0(
      sample,
      "_trimmed_R1.fastq.gz"
    )
  )
  
  read2 <- file.path(
    fastq_dir_abs,
    paste0(
      sample,
      "_trimmed_R2.fastq.gz"
    )
  )
  
  
  # ----------------------------------------------------------
  # Check input
  # ----------------------------------------------------------
  
  if (!file.exists(read1)) {
    
    stop(
      "R1 FASTQ not found:\n",
      read1
    )
    
  }
  
  if (!file.exists(read2)) {
    
    stop(
      "R2 FASTQ not found:\n",
      read2
    )
    
  }
  
  
  cat("\nR1:\n")
  cat(read1, "\n")
  
  cat("\nR2:\n")
  cat(read2, "\n")
  
  
  # ----------------------------------------------------------
  # BAM output prefix
  # ----------------------------------------------------------
  
  output_prefix <- file.path(
    bam_dir_abs,
    paste0(
      sample,
      "_align"
    )
  )
  
  
  # ----------------------------------------------------------
  # Run Subjunc
  # ----------------------------------------------------------
  
  cat("\nStarting RNA-seq alignment...\n")
  cat("This may take a while on a 16 GB laptop.\n\n")
  
  
  alignment <- subjunc(
    
    index = index_base,
    
    readfile1 = read1,
    
    readfile2 = read2,
    
    output_file = paste0(
      output_prefix,
      ".bam"
    ),
    
    input_format = "gzFASTQ",
    
    nthreads = 6,
    
    # Use paired-end fragments
    phredOffset = 33,
    
    # Report uniquely mapped reads in the primary output
    unique = TRUE,
    
    # Detect indels
    indels = 5,
    
    # Allow multiple mapping locations internally;
    # unique filtering controls the primary output.
    nBestLocations = 1,
    
    # Keep junction information
    output_format = "BAM",
  )
  
  
  # ----------------------------------------------------------
  # Print alignment summary
  # ----------------------------------------------------------
  
  cat("\nAlignment complete.\n")
  
  print(alignment)
  
  
  # ----------------------------------------------------------
  # Check BAM
  # ----------------------------------------------------------
  
  bam_file <- paste0(
    output_prefix,
    ".BAM"
  )
  
  # Depending on Rsubread's filename handling, check both
  # possible case conventions.
  
  possible_bams <- c(
    paste0(output_prefix, ".bam"),
    paste0(output_prefix, ".BAM")
  )
  
  existing_bams <- possible_bams[
    file.exists(possible_bams)
  ]
  
  if (length(existing_bams) == 0) {
    
    stop(
      "Alignment completed but no BAM file was found.\n",
      "Expected one of:\n",
      paste(
        possible_bams,
        collapse = "\n"
      )
    )
    
  }
  
  
  bam_file <- existing_bams[1]
  
  
  cat(
    "\nBAM file created:\n",
    bam_file,
    "\n"
  )
  
  cat(
    "BAM size:",
    round(
      file.info(bam_file)$size / 1e9,
      2
    ),
    "GB\n"
  )
  
  
  # ----------------------------------------------------------
  # Save alignment summary
  # ----------------------------------------------------------
  
  saveRDS(
    alignment,
    file.path(
      log_dir,
      paste0(
        sample,
        "_alignment_summary.rds"
      )
    )
  )
  
}


cat(
  "\n========================================\n",
  "All requested alignments complete.\n",
  "========================================\n"
)