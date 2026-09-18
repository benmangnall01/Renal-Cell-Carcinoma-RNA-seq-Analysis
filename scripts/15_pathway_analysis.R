# ============================================================
# 15_pathway_analysis.R
#
# Pathway analysis of pRCC RNA-seq differential expression
#
# Analyses:
#   1. GO Biological Process enrichment - upregulated genes
#   2. GO Biological Process enrichment - downregulated genes
#   3. GO Biological Process GSEA using the full ranked list
#
# Uses org.Hs.eg.db locally.
# ============================================================

library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

setwd("C:/Users/benma/RNA-seq")

# ============================================================
# Directories
# ============================================================

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

# ============================================================
# Load fitted DESeq2 object
# ============================================================

dds <- readRDS(
  file.path(
    processed_dir,
    "dds_fitted.rds"
  )
)

cat("DESeq2 object loaded.\n")

# ============================================================
# Obtain tumour vs normal results
# ============================================================

res <- results(
  dds,
  contrast = c(
    "tissue",
    "Tumour",
    "Normal"
  )
)

res_df <- as.data.frame(res)

res_df$ensembl_id <- rownames(res_df)

# ============================================================
# Clean Ensembl IDs
# ============================================================

res_df <- res_df %>%
  dplyr::mutate(
    ensembl_id_clean = sub(
      "\\..*$",
      "",
      ensembl_id
    )
  )

# ============================================================
# Significant genes
#
# Same threshold as differential_expression.R:
# FDR < 0.05
# absolute log2FC >= 1
# ============================================================

sig_df <- res_df %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )

up_genes <- sig_df %>%
  dplyr::filter(
    log2FoldChange >= 1
  ) %>%
  dplyr::pull(ensembl_id_clean)

down_genes <- sig_df %>%
  dplyr::filter(
    log2FoldChange <= -1
  ) %>%
  dplyr::pull(ensembl_id_clean)

up_genes <- unique(up_genes)
down_genes <- unique(down_genes)

# ============================================================
# Background / universe
# ============================================================

universe <- res_df %>%
  dplyr::filter(
    !is.na(pvalue),
    !is.na(stat)
  ) %>%
  dplyr::pull(ensembl_id_clean)

universe <- unique(universe)

# ============================================================
# Summary
# ============================================================

cat("\n========================================\n")
cat("Pathway analysis summary\n")
cat("========================================\n")

cat(
  "Total genes in DESeq2 results:",
  nrow(res_df),
  "\n"
)

cat(
  "Significant genes:",
  length(unique(sig_df$ensembl_id_clean)),
  "\n"
)

cat(
  "Upregulated genes:",
  length(up_genes),
  "\n"
)

cat(
  "Downregulated genes:",
  length(down_genes),
  "\n"
)

cat(
  "Background genes:",
  length(universe),
  "\n"
)

# ============================================================
# GO enrichment - UPREGULATED
# ============================================================

cat("\nRunning GO enrichment for upregulated genes...\n")

ego_up <- enrichGO(
  gene = up_genes,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ============================================================
# GO enrichment - DOWNREGULATED
# ============================================================

cat("Running GO enrichment for downregulated genes...\n")

ego_down <- enrichGO(
  gene = down_genes,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ============================================================
# Convert results to data frames
# ============================================================

ego_up_df <- as.data.frame(ego_up)
ego_down_df <- as.data.frame(ego_down)

# ============================================================
# Save GO results
# ============================================================

write.csv(
  ego_up_df,
  file.path(
    results_dir,
    "GO_BP_upregulated.csv"
  ),
  row.names = FALSE
)

write.csv(
  ego_down_df,
  file.path(
    results_dir,
    "GO_BP_downregulated.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Print top GO terms - UP
# ============================================================

cat("\n========================================\n")
cat("Top GO Biological Processes - UP\n")
cat("========================================\n")

if (nrow(ego_up_df) > 0) {
  
  top_up <- ego_up_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count
    ) %>%
    head(20)
  
  print(top_up)
  
} else {
  
  cat("No significantly enriched GO terms found.\n")
  
}

# ============================================================
# Print top GO terms - DOWN
# ============================================================

cat("\n========================================\n")
cat("Top GO Biological Processes - DOWN\n")
cat("========================================\n")

if (nrow(ego_down_df) > 0) {
  
  top_down <- ego_down_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count
    ) %>%
    head(20)
  
  print(top_down)
  
} else {
  
  cat("No significantly enriched GO terms found.\n")
  
}

# ============================================================
# Plot GO enrichment - UP
# ============================================================

if (nrow(ego_up_df) > 0) {
  
  p_up <- dotplot(
    ego_up,
    showCategory = 20
  ) +
    ggtitle(
      "GO Biological Process - Upregulated genes"
    )
  
  ggsave(
    file.path(
      figures_dir,
      "06_GO_BP_upregulated.png"
    ),
    p_up,
    width = 10,
    height = 8,
    dpi = 300
  )
  
}

# ============================================================
# Plot GO enrichment - DOWN
# ============================================================

if (nrow(ego_down_df) > 0) {
  
  p_down <- dotplot(
    ego_down,
    showCategory = 20
  ) +
    ggtitle(
      "GO Biological Process - Downregulated genes"
    )
  
  ggsave(
    file.path(
      figures_dir,
      "07_GO_BP_downregulated.png"
    ),
    p_down,
    width = 10,
    height = 8,
    dpi = 300
  )
  
}

# ============================================================
# Prepare ranked list for GSEA
# ============================================================

cat("\nPreparing ranked gene list for GO GSEA...\n")

ranked_genes <- res_df %>%
  dplyr::filter(
    !is.na(stat),
    !is.na(ensembl_id_clean)
  ) %>%
  dplyr::select(
    ensembl_id_clean,
    stat
  )

# ------------------------------------------------------------
# Remove duplicated Ensembl IDs
#
# Keep the entry with the strongest absolute statistic.
# ------------------------------------------------------------

ranked_genes <- ranked_genes %>%
  dplyr::group_by(
    ensembl_id_clean
  ) %>%
  dplyr::slice_max(
    order_by = abs(stat),
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

gene_rank <- ranked_genes$stat

names(gene_rank) <- ranked_genes$ensembl_id_clean

gene_rank <- sort(
  gene_rank,
  decreasing = TRUE
)

cat(
  "Genes in GSEA ranking:",
  length(gene_rank),
  "\n"
)

# ============================================================
# GO GSEA
# ============================================================

cat("\nRunning GO Biological Process GSEA...\n")

gsea_go <- gseGO(
  geneList = gene_rank,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea"
)

# ============================================================
# Save GSEA results
# ============================================================

gsea_go_df <- as.data.frame(gsea_go)

write.csv(
  gsea_go_df,
  file.path(
    results_dir,
    "GSEA_GO_BP.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Print top GSEA pathways
# ============================================================

cat("\n========================================\n")
cat("Top GO Biological Process GSEA pathways\n")
cat("========================================\n")

if (nrow(gsea_go_df) > 0) {
  
  top_gsea <- gsea_go_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(
      ID,
      Description,
      setSize,
      enrichmentScore,
      NES,
      pvalue,
      p.adjust,
      qvalue
    ) %>%
    head(20)
  
} else {
  
  cat("No significantly enriched GO pathways found.\n")
  
}

# ============================================================
# Plot GSEA results
# ============================================================

if (nrow(gsea_go_df) > 0) {
  
  p_gsea <- dotplot(
    gsea_go,
    showCategory = 20
  ) +
    ggtitle(
      "GO Biological Process GSEA"
    )
  
  ggsave(
    file.path(
      figures_dir,
      "08_GSEA_GO_BP.png"
    ),
    p_gsea,
    width = 10,
    height = 8,
    dpi = 300
  )
  
}

# ============================================================
# Finished
# ============================================================

cat("\n========================================\n")
cat("PATHWAY ANALYSIS COMPLETE\n")
cat("========================================\n")

cat("\nResults:\n")

cat(
  "results/GO_BP_upregulated.csv\n"
)

cat(
  "results/GO_BP_downregulated.csv\n"
)

cat(
  "results/GSEA_GO_BP.csv\n"
)

cat("\nFigures:\n")

cat(
  "results/figures/06_GO_BP_upregulated.png\n"
)

cat(
  "results/figures/07_GO_BP_downregulated.png\n"
)

cat(
  "results/figures/08_GSEA_GO_BP.png\n"
)