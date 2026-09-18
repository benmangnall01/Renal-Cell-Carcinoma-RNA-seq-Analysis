# ============================================================
# 07_download_reference.R
#
# Download the GRCh38 primary assembly and the Ensembl 84
# GTF annotation used by the published study.
# ============================================================

setwd("C:/Users/benma/RNA-seq")

library(curl)

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

genome_dir <- file.path(
  "reference",
  "genome"
)

annotation_dir <- file.path(
  "reference",
  "annotation"
)

dir.create(
  genome_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  annotation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Reference URLs
# ------------------------------------------------------------

genome_url <- paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/",
  "Gencode_human/release_28/",
  "GRCh38.primary_assembly.genome.fa.gz"
)

gtf_url <- paste0(
  "https://ftp.ensembl.org/pub/release-84/gtf/",
  "homo_sapiens/",
  "Homo_sapiens.GRCh38.84.gtf.gz"
)


# ------------------------------------------------------------
# Local files
# ------------------------------------------------------------

genome_file <- file.path(
  genome_dir,
  "GRCh38.primary_assembly.genome.fa.gz"
)

gtf_file <- file.path(
  annotation_dir,
  "Homo_sapiens.GRCh38.84.gtf.gz"
)


# ------------------------------------------------------------
# Download helper
# ------------------------------------------------------------

download_if_missing <- function(
    url,
    destination
) {
  
  if (file.exists(destination)) {
    
    cat(
      "\nAlready exists:\n",
      destination,
      "\n"
    )
    
    return(invisible(NULL))
  }
  
  
  cat(
    "\nDownloading:\n",
    url,
    "\n"
  )
  
  cat(
    "To:\n",
    destination,
    "\n\n"
  )
  
  
  curl_download(
    url = url,
    destfile = destination,
    quiet = FALSE
  )
  
  
  if (!file.exists(destination)) {
    
    stop(
      "Download failed:\n",
      destination
    )
  }
  
  
  cat(
    "\nDownloaded successfully.\n"
  )
}


# ------------------------------------------------------------
# Download genome
# ------------------------------------------------------------

download_if_missing(
  genome_url,
  genome_file
)


# ------------------------------------------------------------
# Download annotation
# ------------------------------------------------------------

download_if_missing(
  gtf_url,
  gtf_file
)


# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Reference download complete\n")
cat("========================================\n\n")

cat(
  "Genome:\n",
  genome_file,
  "\n\n"
)

cat(
  "Annotation:\n",
  gtf_file,
  "\n"
)