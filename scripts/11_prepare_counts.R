# Here I switched over to the published gene count matrix's
# Due to not being able to run all 106 4.4gb files on laptop

library(tidyverse)
setwd("C:/Users/benma/RNA-seq")

# -----------------------------
# Directories
# -----------------------------

dataset_dir <- file.path("data", "raw", "GSE180777")
processed_dir <- file.path("data", "processed")

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------
# Load raw count matrix
# -----------------------------

count_file <- file.path(
  dataset_dir,
  "GSE180777_processed_data_files.csv.gz"
)

counts_raw <- read.csv(
  gzfile(count_file),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Raw count matrix dimensions:\n")
print(dim(counts_raw))

# -----------------------------
# Separate gene ID and symbol
# -----------------------------

gene_info <- counts_raw$gene_id

gene_split <- str_split_fixed(
  gene_info,
  "\\|",
  n = 2
)

gene_ids <- gene_split[, 1]
gene_symbols <- gene_split[, 2]

# -----------------------------
# Extract count matrix
# -----------------------------

counts <- counts_raw[, -1]

# Convert to numeric matrix
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

# Add Ensembl IDs as row names
rownames(counts) <- gene_ids

# -----------------------------
# Check sample names
# -----------------------------

metadata_clean <- read.csv(
  file.path(
    processed_dir,
    "sample_metadata_clean.csv"
  ),
  stringsAsFactors = FALSE
)

cat("\nNumber of count columns:",
    ncol(counts), "\n")

cat("Number of metadata samples:",
    nrow(metadata_clean), "\n")

cat("\nSamples in count matrix but not metadata:\n")
print(setdiff(
  colnames(counts),
  metadata_clean$sample
))

cat("\nSamples in metadata but not count matrix:\n")
print(setdiff(
  metadata_clean$sample,
  colnames(counts)
))

# -----------------------------
# Check duplicated Ensembl IDs
# -----------------------------

cat("\nDuplicated Ensembl IDs:\n")
print(sum(duplicated(gene_ids)))

# -----------------------------
# Check missing values
# -----------------------------

cat("\nMissing values:\n")
print(sum(is.na(counts)))

# -----------------------------
# Check that counts are integers
# -----------------------------

cat("\nAny non-integer counts?\n")
print(any(counts %% 1 != 0))

# -----------------------------
# Save gene annotation
# -----------------------------

gene_annotation <- data.frame(
  ensembl_id = gene_ids,
  gene_symbol = gene_symbols
)

write.csv(
  gene_annotation,
  file.path(
    processed_dir,
    "gene_annotation.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# Save cleaned counts
# -----------------------------

saveRDS(
  counts,
  file.path(
    processed_dir,
    "counts_raw.rds"
  )
)

cat("\nCleaned count matrix saved.\n")