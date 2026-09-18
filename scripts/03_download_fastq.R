library(tidyverse)
library(curl)

# ============================================================
# 03_download_fastq.R
#
# Download raw paired-end FASTQ files from ENA
# ============================================================

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

fastq_dir <- file.path(
  "data",
  "raw",
  "fastq",
  "untrimmed"
)

log_dir <- file.path(
  "results",
  "logs"
)

dir.create(
  fastq_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Samples to download
# ------------------------------------------------------------

# Start with ONE sample while testing the pipeline.
samples_to_download <- c("001T")

# ------------------------------------------------------------
# Load GEO metadata
# ------------------------------------------------------------

metadata <- read.csv(
  file.path(
    "data",
    "raw",
    "GSE180777",
    "sample_metadata.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Drop duplicate first column
metadata <- metadata[-1]

# ------------------------------------------------------------
# Extract SRA accessions
# ------------------------------------------------------------

metadata <- metadata %>%
  transmute(
    sample = title,
    geo_accession = geo_accession,
    sra_accession = str_extract(
      relation.1,
      "SRX\\d+"
    )
  )


# ------------------------------------------------------------
# Select samples
# ------------------------------------------------------------

selected <- metadata %>%
  filter(
    sample %in% samples_to_download
  )

if (nrow(selected) != length(samples_to_download)) {
  stop(
    "One or more requested samples could not be found."
  )
}

cat("\nSamples selected:\n")
print(selected)


# ------------------------------------------------------------
# Query ENA
# ------------------------------------------------------------

get_ena_runs <- function(sra_accession) {
  
  ena_url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?",
    "accession=",
    sra_accession,
    "&result=read_run",
    "&fields=run_accession,fastq_ftp,",
    "fastq_md5,fastq_bytes",
    "&format=tsv"
  )
  
  cat(
    "\nQuerying ENA for:",
    sra_accession,
    "\n"
  )
  
  result <- read.delim(
    ena_url,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  result
}


# ------------------------------------------------------------
# Download
# ------------------------------------------------------------

download_manifest <- list()


for (i in seq_len(nrow(selected))) {
  
  sample_name <- selected$sample[i]
  
  sra_accession <- selected$sra_accession[i]
  
  cat("\n========================================\n")
  cat("Sample:", sample_name, "\n")
  cat("SRA:", sra_accession, "\n")
  cat("========================================\n")
  
  
  ena_runs <- get_ena_runs(
    sra_accession
  )
  
  if (nrow(ena_runs) == 0) {
    stop(
      "ENA returned no sequencing runs for ",
      sra_accession
    )
  }
  
  run <- ena_runs[1, ]
  
  cat("\nRun accession:")
  cat(run$run_accession)
  cat("\n")
  
  
  # --------------------------------------------------------
  # Extract FASTQ URLs
  # --------------------------------------------------------
  
  fastq_urls <- unlist(
    strsplit(
      run$fastq_ftp,
      ";"
    )
  )
  
  if (length(fastq_urls) != 2) {
    stop(
      "Expected two paired-end FASTQ files, but found ",
      length(fastq_urls)
    )
  }
  
  
  # --------------------------------------------------------
  # Download R1 + R2
  # --------------------------------------------------------
  
  for (read_number in 1:2) {
    
    # ENA returns the hostname/path without a scheme
    url <- paste0(
      "https://",
      fastq_urls[read_number]
    )
    
    local_file <- file.path(
      fastq_dir,
      paste0(
        sample_name,
        "_R",
        read_number,
        ".fastq.gz"
      )
    )
    
    cat("\nDownloading:\n")
    cat(url, "\n")
    
    cat("To:\n")
    cat(local_file, "\n\n")
    
    
    # ----------------------------------------------------
    # Download using curl
    # ----------------------------------------------------
    
    curl_download(
      url = url,
      destfile = local_file,
      quiet = FALSE
    )
    
    
    # ----------------------------------------------------
    # Check that file exists
    # ----------------------------------------------------
    
    if (!file.exists(local_file)) {
      stop(
        "Download failed: ",
        local_file
      )
    }
    
    downloaded_bytes <- file.info(
      local_file
    )$size
    
    cat(
      "\nDownloaded:",
      round(downloaded_bytes / 1e9, 2),
      "GB\n"
    )
    
    
    # ----------------------------------------------------
    # Save manifest information
    # ----------------------------------------------------
    
    expected_bytes <- as.numeric(
      strsplit(
        run$fastq_bytes,
        ";"
      )[[1]][read_number]
    )
    
    expected_md5 <- strsplit(
      run$fastq_md5,
      ";"
    )[[1]][read_number]
    
    
    download_manifest[[length(download_manifest) + 1]] <- data.frame(
      
      sample = sample_name,
      
      geo_accession =
        selected$geo_accession[i],
      
      sra_accession =
        sra_accession,
      
      run_accession =
        run$run_accession,
      
      read =
        paste0(
          "R",
          read_number
        ),
      
      url =
        url,
      
      local_file =
        local_file,
      
      expected_bytes =
        expected_bytes,
      
      downloaded_bytes =
        downloaded_bytes,
      
      expected_md5 =
        expected_md5,
      
      stringsAsFactors = FALSE
    )
  }
}


# ------------------------------------------------------------
# Save manifest
# ------------------------------------------------------------

download_manifest <- bind_rows(
  download_manifest
)

write.csv(
  download_manifest,
  file.path(
    log_dir,
    "fastq_download_manifest.csv"
  ),
  row.names = FALSE
)


cat("\nDownload complete.\n")