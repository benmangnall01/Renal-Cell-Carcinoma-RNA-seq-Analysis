# ============================================================
# 10_count_genes.R
#
# Assign aligned RNA-seq fragments to genes using
# Rsubread::featureCounts().
#
# Input:
#   001N_align.bam
#
# Output:
#   gene-level raw counts
#
# ============================================================

setwd("C:/Users/benma/RNA-seq")

library(Rsubread)
library(tidyverse)

source("config/config.R")


# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

bam_dir <- file.path(
  "data",
  "processed",
  "bam"
)

annotation_dir <- file.path(
  "reference",
  "annotation"
)

count_dir <- file.path(
  "data",
  "processed",
  "counts"
)

results_dir <- file.path(
  "results"
)


dir.create(
  count_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Absolute paths
# ------------------------------------------------------------

bam_dir_abs <- normalizePath(
  bam_dir,
  winslash = "/",
  mustWork = TRUE
)

gtf_file <- normalizePath(
  file.path(
    annotation_dir,
    "Homo_sapiens.GRCh38.84.gtf.gz"
  ),
  winslash = "/",
  mustWork = TRUE
)

count_dir_abs <- normalizePath(
  count_dir,
  winslash = "/",
  mustWork = TRUE
)


# ------------------------------------------------------------
# Samples
# ------------------------------------------------------------

samples <- TEST_SAMPLES

cat("Samples to count:\n")
print(samples)


# ------------------------------------------------------------
# Find BAM files
# ------------------------------------------------------------

bam_files <- file.path(
  bam_dir_abs,
  paste0(
    samples,
    "_align.bam"
  )
)


# ------------------------------------------------------------
# Check BAM files
# ------------------------------------------------------------

missing_bams <- bam_files[
  !file.exists(bam_files)
]

if (length(missing_bams) > 0) {
  
  stop(
    paste(
      "The following BAM files were not found:",
      paste(
        missing_bams,
        collapse = "\n"
      )
    )
  )
  
}


cat("\nBAM files:\n")
print(bam_files)


# ------------------------------------------------------------
# Check GTF
# ------------------------------------------------------------

cat("\nGTF annotation:\n")
cat(gtf_file, "\n")


# ------------------------------------------------------------
# Run featureCounts
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Running featureCounts\n")
cat("========================================\n\n")

cat(
  "Counting paired-end fragments at gene level.\n"
)

cat(
  "Using GTF gene_id as the gene identifier.\n"
)

cat(
  "Assuming unstranded libraries.\n\n"
)


counts_result <- featureCounts(
  
  # ----------------------------------------------------------
  # Input BAM
  # ----------------------------------------------------------
  
  files = bam_files,
  
  
  # ----------------------------------------------------------
  # Annotation
  # ----------------------------------------------------------
  
  annot.ext = gtf_file,
  
  isGTFAnnotationFile = TRUE,
  
  GTF.featureType = "exon",
  
  GTF.attrType = "gene_id",
  
  
  # ----------------------------------------------------------
  # Gene-level counting
  # ----------------------------------------------------------
  
  useMetaFeatures = TRUE,
  
  
  # ----------------------------------------------------------
  # Paired-end settings
  # ----------------------------------------------------------
  
  isPairedEnd = TRUE,
  
  countReadPairs = TRUE,
  
  # Require both mates to be mapped
  requireBothEndsMapped = TRUE,
  
  
  # Keep standard fragment-length rules
  checkFragLength = TRUE,
  
  minFragLength = 50,
  
  maxFragLength = 600,
  
  
  # ----------------------------------------------------------
  # Multi-mapping
  # ----------------------------------------------------------
  
  countMultiMappingReads = FALSE,
  
  
  # Don't assign one fragment to multiple genes
  allowMultiOverlap = FALSE,
  
  
  # ----------------------------------------------------------
  # Strandness
  # ----------------------------------------------------------
  
  # The public metadata available to us does not establish
  # that the library is stranded.
  #
  # 0 = unstranded
  # 1 = stranded
  # 2 = reversely stranded
  
  strandSpecific = 0,
  
  
  # ----------------------------------------------------------
  # Performance
  # ----------------------------------------------------------
  
  nthreads = 6,
  
  verbose = TRUE
  
)


# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("featureCounts complete\n")
cat("========================================\n\n")


cat("Number of genes:\n")
print(
  nrow(
    counts_result$counts
  )
)


cat("\nNumber of samples:\n")
print(
  ncol(
    counts_result$counts
  )
)


# ------------------------------------------------------------
# Assignment statistics
# ------------------------------------------------------------

cat("\nAssignment statistics:\n")

print(
  counts_result$stat
)


# ------------------------------------------------------------
# Save complete featureCounts object
# ------------------------------------------------------------

saveRDS(
  
  counts_result,
  
  file.path(
    count_dir_abs,
    "featureCounts_result.rds"
  )
  
)


# ------------------------------------------------------------
# Extract count matrix
# ------------------------------------------------------------

gene_counts <- counts_result$counts


# ------------------------------------------------------------
# Add sample names
# ------------------------------------------------------------

colnames(gene_counts) <- samples


# ------------------------------------------------------------
# Save raw count matrix
# ------------------------------------------------------------

saveRDS(
  
  gene_counts,
  
  file.path(
    count_dir_abs,
    "001N_gene_counts.rds"
  )
  
)


# ------------------------------------------------------------
# Save as CSV
# ------------------------------------------------------------

gene_counts_df <- as.data.frame(
  gene_counts
)

gene_counts_df$gene_id <- rownames(
  gene_counts
)

gene_counts_df <- gene_counts_df %>%
  select(
    gene_id,
    everything()
  )


write.csv(
  
  gene_counts_df,
  
  file.path(
    count_dir_abs,
    "001N_gene_counts.csv"
  ),
  
  row.names = FALSE
)


# ------------------------------------------------------------
# Save assignment statistics
# ------------------------------------------------------------

write.csv(
  
  counts_result$stat,
  
  file.path(
    results_dir,
    "featureCounts_assignment_stats_001N.csv"
  ),
  
  row.names = FALSE
)


# ------------------------------------------------------------
# Show first genes
# ------------------------------------------------------------

cat("\nFirst ten genes:\n")

print(
  head(
    gene_counts_df,
    10
  )
)

# ------------------------------------------------------------
# Total assigned fragments
# ------------------------------------------------------------

cat("\nTotal assigned fragments:\n")

print(sum(gene_counts))

cat("Gene counting complete.\n")
