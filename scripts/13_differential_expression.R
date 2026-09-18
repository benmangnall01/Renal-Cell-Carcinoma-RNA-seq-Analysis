library(DESeq2)
library(tidyverse)
setwd("C:/Users/benma/RNA-seq")

# -----------------------------
# Directories
# -----------------------------

processed_dir <- file.path("data", "processed")
results_dir <- file.path("results")

# -----------------------------
# Load filtered counts
# -----------------------------

counts_filtered <- readRDS(file.path(processed_dir, "counts_filtered.rds"))

# -----------------------------
# Load metadata
# -----------------------------

metadata <- read.csv(file.path(processed_dir, "sample_metadata_clean.csv"), stringsAsFactors = FALSE)

# Put metadata in exactly the same order as count columns
metadata <- metadata[match(colnames(counts_filtered), metadata$sample),]

rownames(metadata) <- metadata$sample

metadata$tissue <- factor(metadata$tissue, levels = c("Normal", "Tumour"))

metadata$patient <- factor(metadata$patient)

# Safety check
stopifnot(all(rownames(metadata) == colnames(counts_filtered)))

# -----------------------------
# Create DESeq2 object
# -----------------------------

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ patient + tissue)

# -----------------------------
# Run DESeq2
# -----------------------------

dds <- DESeq(dds)

# Save complete DESeq2 object
saveRDS(dds, file.path(processed_dir, "dds_fitted.rds"))

# -----------------------------
# Inspect model coefficients
# -----------------------------

cat("\nDESeq2 coefficients:\n")
print(resultsNames(dds))

# -----------------------------
# Tumour vs Normal
# -----------------------------

res <- results(dds, contrast = c("tissue", "Tumour", "Normal"))

# Convert to data frame
res_df <- as.data.frame(res)

# Add Ensembl ID
res_df$ensembl_id <- rownames(res_df)

# -----------------------------
# Add gene symbols
# -----------------------------

gene_annotation <- read.csv(file.path(processed_dir, "gene_annotation.csv"), stringsAsFactors = FALSE)

res_df <- res_df %>%
  left_join(gene_annotation, by = "ensembl_id")

# -----------------------------
# Order by adjusted p-value
# -----------------------------

res_df <- res_df %>%
  arrange(padj)

# -----------------------------
# Save complete results
# -----------------------------

write.csv(res_df, file.path(results_dir, "DESeq2_tumour_vs_normal_all.csv"), row.names = FALSE)

# -----------------------------
# Significant genes
# -----------------------------

# Standard starting threshold:
# FDR < 0.05 and absolute log2FC >= 1
# This corresponds to at least a 2-fold change.
significant <- res_df %>%
  filter(!is.na(padj), padj < 0.05,abs(log2FoldChange) >= 1)

write.csv(significant, file.path(results_dir, "DESeq2_tumour_vs_normal_significant.csv"), row.names = FALSE)

# -----------------------------
# Summary
# -----------------------------

cat("\nTotal genes tested:", nrow(res_df), "\n")

cat("Genes with adjusted p-value < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")

cat("Genes with FDR < 0.05 and |log2FC| >= 1:", nrow(significant), "\n")

cat("\nTop 20 genes:\n")

print(
  significant %>%
    select(ensembl_id, gene_symbol, baseMean, log2FoldChange, pvalue, padj) %>%
    head(20))