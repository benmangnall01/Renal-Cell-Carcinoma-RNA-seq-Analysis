# ============================================================
# 06_fastqc_trimmed.R
#
# Run FastQC on trimmed paired-end FASTQ files.
# ============================================================

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
  "trimmed"
)

output_dir <- file.path(
  "qc",
  "fastqc_trimmed"
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
# Absolute paths
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
      FASTQC_EXE
    )
  )
  
}


# ------------------------------------------------------------
# Find trimmed FASTQ files
# ------------------------------------------------------------

fastq_files <- list.files(
  input_dir_abs,
  pattern = "\\.fastq\\.gz$",
  full.names = TRUE
)


if (length(fastq_files) == 0) {
  
  stop(
    "No trimmed FASTQ files were found in:\n",
    input_dir_abs
  )
  
}


cat("Trimmed FASTQ files found:\n")

print(
  fastq_files
)


# ------------------------------------------------------------
# Ensure absolute paths
# ------------------------------------------------------------

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
  
  cat(
    "Running FastQC on:\n",
    basename(fastq_file),
    "\n"
  )
  
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
        
        "-classpath",
        ".;./sam-1.103.jar;./jbzip2-0.9.jar",
        
        "uk.ac.babraham.FastQC.FastQCApplication",
        
        fastq_file
        
      ),
      
      stdout = TRUE,
      
      stderr = TRUE
      
    )
    
  )
  
  
  cat(
    paste(
      result,
      collapse = "\n"
    )
  )
  
  cat("\n")
  
  
  # ----------------------------------------------------------
  # Check exit status
  # ----------------------------------------------------------
  
  status <- attr(
    result,
    "status"
  )
  
  
  if (!is.null(status) && status != 0) {
    
    stop(
      paste0(
        "FastQC failed for ",
        basename(fastq_file),
        "\nExit status: ",
        status
      )
    )
    
  }
  
}


cat(
  "\nFastQC analysis of trimmed reads complete.\n"
)

cat(
  "\nReports saved to:\n",
  output_dir_abs,
  "\n"
)