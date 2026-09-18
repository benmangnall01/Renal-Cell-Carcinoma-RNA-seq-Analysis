setwd("C:/Users/benma/RNA-seq")
library(tidyverse)
library(withr)

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
  "qc",
  "fastqc_raw"
)

log_dir <- file.path(
  "results",
  "logs"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Convert paths to absolute paths BEFORE changing directory
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

# ------------------------------------------------------------
# Check FastQC
# ------------------------------------------------------------

if (!file.exists(FASTQC_EXE)) {
  
  stop(
    paste0(
      "FastQC executable was not found at:\n",
      FASTQC_EXE,
      "\n"
    )
  )
}

cat("FastQC executable:\n")
cat(FASTQC_EXE, "\n\n")

# ------------------------------------------------------------
# Test FastQC
# ------------------------------------------------------------

cat("Testing FastQC installation...\n")

fastqc_test <- withr::with_dir(
  FASTQC_DIR,
  system2(
    command = FASTQC_EXE,
    args = "--version",
    stdout = TRUE,
    stderr = TRUE
  )
)

print(fastqc_test)

# ------------------------------------------------------------
# Find FASTQ files
# ------------------------------------------------------------

fastq_files <- list.files(
  input_dir_abs,
  pattern = "\\.fastq\\.gz$",
  full.names = TRUE
)

if (length(fastq_files) == 0) {
  
  stop(
    "No FASTQ files found in:\n",
    input_dir_abs
  )
}

cat("\nFASTQ files found:\n")
print(fastq_files)

# Convert every FASTQ path to absolute path explicitly
fastq_files <- normalizePath(
  fastq_files,
  winslash = "/",
  mustWork = TRUE
)

# ------------------------------------------------------------
# Run FastQC
# ------------------------------------------------------------

for (fastq_file in fastq_files) {
  
  cat("\n========================================\n")
  cat("Running FastQC on:\n")
  cat(basename(fastq_file), "\n")
  cat("========================================\n")
  
  result <- withr::with_dir(
    FASTQC_DIR,
    
    system2(
      command = FASTQC_JAVA,
      
      args = c(
        "-Xmx250m",

        paste0(
          "-Dfastqc.output_dir=",
          output_dir_abs
        ),
        
        # FastQC dependencies
        "-classpath",
        ".;./sam-1.103.jar;./jbzip2-0.9.jar",
        
        # Main FastQC Java class
        "uk.ac.babraham.FastQC.FastQCApplication",
        
        # Input FASTQ
        fastq_file
      ),
      
      stdout = TRUE,
      stderr = TRUE
    )
  )
  
  cat(paste(result, collapse = "\n"))
  
  cat("\n")
  
  
  # --------------------------------------------------------
  # Check exit status
  # --------------------------------------------------------
  
  status <- attr(result, "status")
  
  if (!is.null(status) && status != 0) {
    stop(paste0("FastQC failed for ", basename(fastq_file), "\nExit status: ", status))
  }
}

cat("\nFastQC processing complete.\n")

cat(
  "Reports saved to:\n",
  output_dir_abs,
  "\n"
)