library(GEOquery)

# -----------------------------
# Project directories
# -----------------------------

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Dataset
# -----------------------------

geo_id <- "GSE180777"

# Download supplementary files from GEO
getGEOSuppFiles(
  geo_id,
  baseDir = raw_dir,
  makeDirectory = TRUE
)

cat("Download complete.\n")

# -----------------------------
# Download sample metadata
# -----------------------------

# Download GEO metadata
geo_data <- getGEO(geo_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

# There is one Series Matrix object for this study
eset <- geo_data[[1]]

metadata <- pData(eset)

write.csv(metadata, file.path(raw_dir, geo_id, "sample_metadata.csv"), row.names = TRUE)

cat("Metadata saved.\n")

# -----------------------------
# Inspect the data
# -----------------------------

dataset_dir <- file.path(raw_dir, geo_id)

print(list.files(dataset_dir))

metadata <- read.csv(
  file.path(dataset_dir, "sample_metadata.csv"),
  check.names = FALSE
)

cat("\nMetadata dimensions:\n")
print(dim(metadata))

cat("\nMetadata columns:\n")
print(colnames(metadata))

cat("\nFirst few rows:\n")
print(head(metadata))

count_file <- file.path(
  dataset_dir,
  "GSE180777_processed_data_files.csv.gz"
)

counts_preview <- read.csv(
  gzfile(count_file),
  nrows = 5,
  check.names = FALSE
)

cat("\nCount-file dimensions (preview only):\n")
print(dim(counts_preview))

cat("\nColumn names:\n")
print(colnames(counts_preview))

cat("\nPreview:\n")
print(counts_preview)