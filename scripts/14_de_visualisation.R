library(DESeq2)
library(tidyverse)

setwd("C:/Users/benma/RNA-seq")

processed_dir <- file.path("data", "processed")
results_dir <- file.path("results")
figures_dir <- file.path(results_dir, "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load fitted DESeq2 object
# ------------------------------------------------------------

dds <- readRDS(
  file.path(
    processed_dir,
    "dds_fitted.rds"
  )
)

# ------------------------------------------------------------
# Obtain tumour vs normal results
# ------------------------------------------------------------

res <- results(
  dds,
  contrast = c("tissue", "Tumour", "Normal")
)

res_df <- as.data.frame(res)

res_df$ensembl_id <- rownames(res_df)

gene_annotation <- read.csv(
  file.path(
    processed_dir,
    "gene_annotation.csv"
  ),
  stringsAsFactors = FALSE
)

res_df <- res_df %>%
  left_join(
    gene_annotation,
    by = "ensembl_id"
  )

# ------------------------------------------------------------
# MA plot
# ------------------------------------------------------------

png(
  file.path(figures_dir, "04_MA_plot.png"),
  width = 2400,
  height = 1800,
  res = 300
)

plotMA(
  res,
  alpha = 0.05,
  ylim = c(-5, 5),
  main = "Tumour vs Normal RNA-seq"
)

dev.off()

# ------------------------------------------------------------
# Prepare volcano plot
# ------------------------------------------------------------

volcano_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = case_when(
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange >= 1 ~ "Up",
      
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange <= -1 ~ "Down",
      
      TRUE ~ "Not significant"
    )
  )

# Avoid infinite values
volcano_df$neg_log10_padj[
  is.infinite(volcano_df$neg_log10_padj)
] <- NA

# ------------------------------------------------------------
# Volcano plot
# ------------------------------------------------------------

p <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj
  )
) +
  geom_point(
    aes(colour = significance),
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential expression: Tumour vs Normal",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    colour = "Category"
  ) +
  theme_minimal()

ggsave(
  file.path(
    figures_dir,
    "05_volcano_plot.png"
  ),
  p,
  width = 9,
  height = 7,
  dpi = 300
)

cat("\nDE visualisation complete.\n")