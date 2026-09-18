library(tidyverse)

# -----------------------------
# Paths
# -----------------------------

dataset_dir <- file.path("data", "raw", "GSE180777")

metadata_file <- file.path(
  dataset_dir,
  "sample_metadata.csv"
)

# -----------------------------
# Load GEO metadata
# -----------------------------

metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# The first column is the CSV row-name column
metadata <- metadata[, -1]

# -----------------------------
# Create our own sample labels
# -----------------------------

metadata_clean <- metadata %>%
  transmute(
    sample = title,
    geo_accession = geo_accession
  ) %>%
  mutate(
    patient = str_extract(sample, "^\\d+"),
    tissue = case_when(
      str_detect(sample, "N$") ~ "Normal",
      str_detect(sample, "T$") ~ "Tumour",
      TRUE ~ NA_character_
    )
  )

# Make tissue a factor and explicitly set the reference
metadata_clean$tissue <- factor(
  metadata_clean$tissue,
  levels = c("Normal", "Tumour")
)

# -----------------------------
# Validation
# -----------------------------

cat("Number of samples:", nrow(metadata_clean), "\n")

cat("\nTissue counts:\n")
print(table(metadata_clean$tissue))

cat("\nNumber of unique patients:",
    length(unique(metadata_clean$patient)), "\n")

cat("\nMissing values:\n")
print(colSums(is.na(metadata_clean)))

cat("\nFirst 10 samples:\n")
print(head(metadata_clean, 10))

# -----------------------------
# Check each patient has N + T
# -----------------------------

pair_check <- metadata_clean %>%
  count(patient, tissue) %>%
  pivot_wider(
    names_from = tissue,
    values_from = n,
    values_fill = 0
  )

cat("\nPairing check:\n")
print(pair_check)

# -----------------------------
# Save
# -----------------------------

write.csv(metadata_clean, file.path("data", "processed", "sample_metadata_clean.csv"), row.names = FALSE)