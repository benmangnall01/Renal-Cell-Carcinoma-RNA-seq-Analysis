library(DESeq2)
library(tidyverse)
setwd("C:/Users/benma/RNA-seq")

# -----------------------------
# Directories
# -----------------------------

processed_dir <- file.path("data", "processed")
results_dir <- file.path("results")
figures_dir <- file.path(results_dir, "figures")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------
# Load data
# -----------------------------

counts <- readRDS(
  file.path(
    processed_dir,
    "counts_raw.rds"
  )
)

metadata <- read.csv(
  file.path(
    processed_dir,
    "sample_metadata_clean.csv"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# Make absolutely sure that
# metadata and count matrix
# have identical sample order
# -----------------------------

if (!all(colnames(counts) %in% metadata$sample)) {
  stop("Some count-matrix samples are missing from metadata.")
}

metadata <- metadata[
  match(colnames(counts), metadata$sample),
]

if (!all(metadata$sample == colnames(counts))) {
  stop("Sample ordering does not match.")
}

rownames(metadata) <- metadata$sample

metadata$tissue <- factor(
  metadata$tissue,
  levels = c("Normal", "Tumour")
)

metadata$patient <- factor(
  metadata$patient
)

cat("Samples:", ncol(counts), "\n")
cat("Genes:", nrow(counts), "\n")

# -----------------------------
# 1. Library size
# -----------------------------

library_size <- colSums(counts)

library_qc <- data.frame(
  sample = names(library_size),
  library_size = as.numeric(library_size),
  tissue = metadata[names(library_size), "tissue"],
  patient = metadata[names(library_size), "patient"]
)

cat("\nLibrary size summary:\n")
print(summary(library_qc$library_size))

write.csv(
  library_qc,
  file.path(results_dir, "library_size_qc.csv"),
  row.names = FALSE
)

# Plot library size

p1 <- ggplot(
  library_qc,
  aes(
    x = reorder(sample, library_size),
    y = library_size,
    colour = tissue
  )
) +
  geom_point(size = 2) +
  labs(
    title = "RNA-seq library size",
    x = "Sample",
    y = "Total number of reads"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  file.path(figures_dir, "01_library_size.png"),
  p1,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------
# 2. Number of detected genes
# -----------------------------

detected_genes <- colSums(counts > 0)

detected_qc <- data.frame(
  sample = names(detected_genes),
  detected_genes = as.numeric(detected_genes),
  tissue = metadata[names(detected_genes), "tissue"],
  patient = metadata[names(detected_genes), "patient"]
)

cat("\nDetected genes per sample:\n")
print(summary(detected_qc$detected_genes))

write.csv(
  detected_qc,
  file.path(results_dir, "detected_genes_qc.csv"),
  row.names = FALSE
)

p2 <- ggplot(
  detected_qc,
  aes(
    x = reorder(sample, detected_genes),
    y = detected_genes,
    colour = tissue
  )
) +
  geom_point(size = 2) +
  labs(
    title = "Number of detected genes per sample",
    x = "Sample",
    y = "Genes with count > 0"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  file.path(figures_dir, "02_detected_genes.png"),
  p2,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------
# 3. Filter very low-count genes
# -----------------------------
#
# Keep genes with at least 10 reads
# in at least 10 samples.
#
# This removes genes that contain
# almost no information for this
# experiment.

keep <- rowSums(counts >= 10) >= 10

counts_filtered <- counts[keep, ]

cat("\nGenes before filtering:",
    nrow(counts), "\n")

cat("Genes after filtering:",
    nrow(counts_filtered), "\n")

cat("Genes removed:",
    nrow(counts) - nrow(counts_filtered), "\n")

saveRDS(
  counts_filtered,
  file.path(
    processed_dir,
    "counts_filtered.rds"
  )
)

# -----------------------------
# 4. Create DESeq2 object
# -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ patient + tissue
)

# -----------------------------
# 5. Variance stabilising
#    transformation
# -----------------------------

vsd <- vst(
  dds,
  blind = TRUE
)

saveRDS(
  vsd,
  file.path(
    processed_dir,
    "vsd.rds"
  )
)

# -----------------------------
# 6. PCA
# -----------------------------

# Use the 500 most variable genes.
vsd_matrix <- assay(vsd)

gene_variances <- apply(
  vsd_matrix,
  1,
  var
)

top_n <- min(
  500,
  length(gene_variances)
)

top_genes <- names(
  sort(
    gene_variances,
    decreasing = TRUE
  )
)[1:top_n]

pca <- prcomp(
  t(vsd_matrix[top_genes, ]),
  scale. = FALSE
)

variance_explained <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_df <- left_join(
  pca_df,
  metadata,
  by = "sample"
)

write.csv(
  pca_df,
  file.path(
    results_dir,
    "pca_coordinates.csv"
  ),
  row.names = FALSE
)

cat("\nVariance explained:\n")

cat(
  "PC1:",
  round(100 * variance_explained[1], 2),
  "%\n"
)

cat(
  "PC2:",
  round(100 * variance_explained[2], 2),
  "%\n"
)

# PCA plot

p3 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = tissue
  )
) +
  geom_line(
    aes(group = patient),
    alpha = 0.25,
    colour = "grey50"
  ) +
  geom_point(size = 3) +
  labs(
    title = "PCA of RNA-seq samples",
    x = paste0(
      "PC1 (",
      round(100 * variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(100 * variance_explained[2], 1),
      "%)"
    )
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, "03_PCA.png"),
  p3,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nQC analysis complete.\n")