# ============================================================
# 05_trim_reads.R
# ============================================================

# ------------------------------------------------------------
# Project directory
# ------------------------------------------------------------

setwd("C:/Users/benma/RNA-seq")

# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------

library(Rfastp)
library(tidyverse)

# ------------------------------------------------------------
# Project configuration
# ------------------------------------------------------------

source("config/config.R")

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

input_dir <- file.path(
  "data",
  "raw",
  "fastq",
  "untrimmed"
)

output_dir <- file.path(
  "data",
  "raw",
  "fastq",
  "trimmed"
)

report_dir <- file.path(
  "qc",
  "rfastp"
)

results_dir <- file.path(
  "results"
)

log_dir <- file.path(
  results_dir,
  "logs"
)


# ------------------------------------------------------------
# Create directories
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  report_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Convert directory paths to absolute paths
# ------------------------------------------------------------

input_dir_abs <- normalizePath(
  input_dir,
  winslash = "/",
  mustWork = TRUE
)

output_dir_abs <- normalizePath(
  output_dir,
  winslash = "/",
  mustWork = TRUE
)

report_dir_abs <- normalizePath(
  report_dir,
  winslash = "/",
  mustWork = TRUE
)

results_dir_abs <- normalizePath(
  results_dir,
  winslash = "/",
  mustWork = TRUE
)

log_dir_abs <- normalizePath(
  log_dir,
  winslash = "/",
  mustWork = TRUE
)


# ------------------------------------------------------------
# Check Rfastp version
# ------------------------------------------------------------

cat("========================================\n")
cat("Rfastp RNA-seq preprocessing\n")
cat("========================================\n\n")

cat(
  "Rfastp version: ",
  as.character(
    packageVersion("Rfastp")
  ),
  "\n\n"
)

# ------------------------------------------------------------
# Samples
# ------------------------------------------------------------

samples <- TEST_SAMPLES

cat("Samples to process:\n")
print(samples)

# ------------------------------------------------------------
# Storage for results
# ------------------------------------------------------------

manifest <- list()

qc_summaries <- list()

trim_summaries <- list()

# ------------------------------------------------------------
# Process samples
# ------------------------------------------------------------

for (sample in samples) {
  
  cat("\n========================================\n")
  cat("Processing sample:", sample, "\n")
  cat("========================================\n")
  
  # ----------------------------------------------------------
  # Input FASTQ files
  # ----------------------------------------------------------
  
  input_r1 <- file.path(
    input_dir_abs,
    paste0(
      sample,
      "_R1.fastq.gz"
    )
  )
  
  input_r2 <- file.path(
    input_dir_abs,
    paste0(
      sample,
      "_R2.fastq.gz"
    )
  )
  
  
  # ----------------------------------------------------------
  # Check input files
  # ----------------------------------------------------------
  
  cat("\nChecking input files...\n")
  
  cat(
    "R1:\n",
    input_r1,
    "\n"
  )
  
  cat(
    "R2:\n",
    input_r2,
    "\n"
  )
  
  
  if (!file.exists(input_r1)) {
    
    stop(
      "R1 does not exist:\n",
      input_r1
    )
    
  }
  
  
  if (!file.exists(input_r2)) {
    
    stop(
      "R2 does not exist:\n",
      input_r2
    )
    
  }
  
  
  cat("Input files found.\n")
  
  
  # ----------------------------------------------------------
  # Output prefix
  # ----------------------------------------------------------
  #
  # Rfastp uses outputFastq as a PREFIX.
  #
  # For example:
  #
  # outputFastq =
  #   data/raw/fastq/trimmed/001N_trimmed
  #
  # produces:
  #
  # 001N_trimmed_R1.fastq.gz
  # 001N_trimmed_R2.fastq.gz
  #
  # ----------------------------------------------------------
  
  output_prefix <- file.path(
    output_dir_abs,
    paste0(
      sample,
      "_trimmed"
    )
  )
  
  
  output_r1 <- paste0(
    output_prefix,
    "_R1.fastq.gz"
  )
  
  output_r2 <- paste0(
    output_prefix,
    "_R2.fastq.gz"
  )
  
  
  # ----------------------------------------------------------
  # Rfastp report location
  #
  # The report object is returned directly by rfastp().
  #
  # We save that object as an RDS after the run.
  # ----------------------------------------------------------
  
  report_rds <- file.path(
    report_dir_abs,
    paste0(
      sample,
      "_rfastp_report.rds"
    )
  )
  
  
  # ----------------------------------------------------------
  # Run Rfastp
  # ----------------------------------------------------------
  #
  # Our FastQC results showed:
  #
  #   - good per-base quality
  #   - no meaningful adapter-content signal
  #   - no N-content problem
  #
  # Therefore we do NOT forcibly remove bases from the
  # front or tail of the reads.
  #
  # We allow adapter detection/trimming and conservative
  # quality/length filtering.
  #
  # ----------------------------------------------------------
  
  cat("\nRunning Rfastp...\n\n")
  
  
  rfastp_report <- rfastp(
    
    # paired-end input
    read1 = input_r1,
    
    read2 = input_r2,
    
    
    # output prefix
    outputFastq = output_prefix,
    
    adapterTrimming = TRUE,
    
    adapterSequenceRead1 = "auto",
    
    adapterSequenceRead2 = "auto",
    
    trimFrontRead1 = 0,
    
    trimTailRead1 = 0,
    
    trimFrontRead2 = 0,
    
    trimTailRead2 = 0,
    
    cutLowQualFront = FALSE,
    
    cutLowQualTail = FALSE,
    
    cutSlideWindowRight = FALSE,
    
    qualityFiltering = TRUE,
    
    qualityFilterPhred = 15,
    
    qualityFilterPercent = 40,
    
    maxNfilter = 5, # Retain the default maximum of 5 Ns.
    
    lengthFiltering = TRUE,
    
    minReadLength = 30,
    
    umi = FALSE,
    
    overrepresentationAnalysis = FALSE,
    
    thread = 8,
    
    verbose = FALSE # Suppress extremely verbose output
  )
  
  # ----------------------------------------------------------
  # Save the Rfastp report object
  # ----------------------------------------------------------
  
  saveRDS(rfastp_report, report_rds)
  
  # ----------------------------------------------------------
  # Generate summaries
  # ----------------------------------------------------------
  
  qc <- qcSummary(
    rfastp_report
  )
  
  trim <- trimSummary(
    rfastp_report
  )
  
  
  # Store results
  
  qc_summaries[[sample]] <- qc
  
  trim_summaries[[sample]] <- trim
  
  # ----------------------------------------------------------
  # Print summaries
  # ----------------------------------------------------------
  
  cat("\nQC summary:\n")
  print(qc)
  
  cat("\nTrim summary:\n")
  print(trim)
  
  
  # ----------------------------------------------------------
  # Check expected FASTQ output
  # ----------------------------------------------------------
  
  expected_outputs <- c(
    output_r1,
    output_r2
  )
  
  
  missing_outputs <- expected_outputs[
    !file.exists(expected_outputs)
  ]
  
  
  if (length(missing_outputs) > 0) {
    
    stop(
      paste(
        "\nRfastp did not create the expected FASTQ files:",
        paste(
          missing_outputs,
          collapse = "\n"
        )
      )
    )
    
  }
  
  
  cat(
    "\nTrimmed FASTQ files created successfully.\n"
  )
  
  
  # ----------------------------------------------------------
  # Record manifest entry
  # ----------------------------------------------------------
  
  manifest[[length(manifest) + 1]] <- data.frame(
    
    sample = sample,
    
    input_r1 = input_r1,
    
    input_r2 = input_r2,
    
    output_r1 = output_r1,
    
    output_r2 = output_r2,
    
    report_rds = report_rds,
    
    stringsAsFactors = FALSE
    
  )
  
  
  cat(
    "\nCompleted sample:",
    sample,
    "\n"
  )
  
}


# ------------------------------------------------------------
# Combine manifest
# ------------------------------------------------------------

manifest <- bind_rows(
  manifest
)


# ------------------------------------------------------------
# Save manifest
# ------------------------------------------------------------

write.csv(
  
  manifest,
  
  file.path(
    log_dir_abs,
    "rfastp_manifest.csv"
  ),
  
  row.names = FALSE
)


# ------------------------------------------------------------
# Combine QC summaries
# ------------------------------------------------------------

qc_combined <- bind_rows(
  qc_summaries,
  .id = "sample"
)

trim_combined <- bind_rows(
  trim_summaries,
  .id = "sample"
)


# ------------------------------------------------------------
# Save summaries
# ------------------------------------------------------------

write.csv(
  
  qc_combined,
  
  file.path(
    results_dir_abs,
    "rfastp_qc_summary.csv"
  ),
  
  row.names = FALSE
  
)


write.csv(
  
  trim_combined,
  
  file.path(
    results_dir_abs,
    "rfastp_trim_summary.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# Final message
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Rfastp processing complete.\n")
cat("========================================\n\n")

cat(
  "Trimmed FASTQ directory:\n",
  output_dir_abs,
  "\n\n"
)

cat(
  "Rfastp reports:\n",
  report_dir_abs,
  "\n\n"
)

cat(
  "QC summary:\n",
  file.path(
    results_dir_abs,
    "rfastp_qc_summary.csv"
  ),
  "\n\n"
)

cat(
  "Trim summary:\n",
  file.path(
    results_dir_abs,
    "rfastp_trim_summary.csv"
  ),
  "\n"
)