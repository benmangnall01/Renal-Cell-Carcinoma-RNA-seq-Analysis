# ============================================================
# 08_build_reference_index.R
#
# Build a gapped Rsubread index for GRCh38.
# ============================================================

setwd("C:/Users/benma/RNA-seq")

library(Rsubread)

# ------------------------------------------------------------
# Files
# ------------------------------------------------------------

genome_file <- normalizePath(
  "reference/genome/GRCh38.primary_assembly.genome.fa.gz",
  winslash = "/",
  mustWork = TRUE
)

index_dir <- normalizePath(
  "reference/genome",
  winslash = "/",
  mustWork = TRUE
)

index_base <- file.path(
  index_dir,
  "GRCh38_rsubread"
)


# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

cat("========================================\n")
cat("Rsubread genome index\n")
cat("========================================\n\n")

cat("Rsubread version:\n")
print(packageVersion("Rsubread"))

cat("\nReference genome:\n")
cat(genome_file, "\n")

cat("\nIndex base:\n")
cat(index_base, "\n")


# ------------------------------------------------------------
# Build index
# ------------------------------------------------------------

cat("\nBuilding gapped index...\n")
cat("This is a one-off operation and may take some time.\n\n")


buildindex(
  basename = index_base,
  
  reference = genome_file,
  
  # Gapped index is much more memory-efficient
  # than a full index.
  gappedIndex = TRUE,
  
  # Keep the index as one block where possible.
  # Mapping memory is controlled separately.
  indexSplit = FALSE,
  
  # Memory available to the index during mapping.
  #
  # 6000 MB is deliberately conservative for a
  # 16 GB laptop.
  memory = 6000,
  
  # Default repeat threshold
  TH_subread = 100
)


# ------------------------------------------------------------
# Check generated index
# ------------------------------------------------------------

index_files <- list.files(
  index_dir,
  pattern = "^GRCh38_rsubread",
  full.names = TRUE
)

cat("\nGenerated index files:\n")
print(index_files)


if (length(index_files) == 0) {
  
  stop(
    "No Rsubread index files were generated."
  )
  
}


cat(
  "\nReference index successfully created.\n"
)