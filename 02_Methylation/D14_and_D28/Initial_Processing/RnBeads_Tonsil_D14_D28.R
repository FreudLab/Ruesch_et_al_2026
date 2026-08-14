# Get the temporary data directory from command line arguments
args <- commandArgs(trailingOnly = TRUE)
data.dir <- args[1]

# Load necessary libraries
library(RnBeads)
library(RnBeads.hg38)
library(grid)

# Check package version (ensure compatibility)
packageVersion("RnBeads")  # Should match your needs (2.16.0 for EPIC v1)

# Define directories
idat.dir <- file.path(data.dir, "idats")
sample.annotation <- file.path(data.dir, "sample_sheet_tonsil_NKDI_and_in_vitro.csv")
analysis.dir <- file.path(data.dir, "RnBeads_Output")
report.dir <- file.path(analysis.dir, "reports")

# ===== CRITICAL FIXES =====
# 1. Verify sample sheet and IDAT files
sample_info <- read.csv(sample.annotation)
print("Sample sheet head:")
print(head(sample_info))

# Check for missing IDAT files
expected_idat <- paste0(sample_info$Sentrix_ID, "_", sample_info$Sentrix_Position, "_Grn.idat")
missing_idat <- expected_idat[!file.exists(file.path(idat.dir, expected_idat))]
if (length(missing_idat) > 0) {
  stop("Missing IDAT files: ", paste(missing_idat, collapse = ", "))
} else {
  print("All IDAT files are present.")
}

# ===== RECOMMENDED CHANGES =====
# Use Sentrix_Position as the primary identifier (more reliable)
rnb.options(
  identifiers.column = "Sample_ID",
  differential = TRUE,
  exploratory = FALSE,
  region.types = NULL,
  normalization.background.method = "none",
  normalization.method = "bmiq",
  filtering.context.removal = NULL,
  filtering.snp = "no",
  filtering.cross.reactive = FALSE,
  filtering.sex.chromosomes.removal = FALSE,
  filtering.missing.value.quantile = 1,
  filtering.coverage.threshold = 0,
  filtering.greedycut = FALSE,
  assembly = "hg38",
  disk.dump.big.matrices = FALSE,
  export.to.csv = TRUE,
  export.to.bed = TRUE,
  export.to.trackhub = NULL,
  import.bed.style = "BisSNP"
)

# Enable parallel processing
logger.start(fname = NA)
num.cores <- 20
parallel.setup(num.cores)
if (!parallel.isEnabled()) {
  warning("Parallel processing failed. Running in serial mode.")
}

# ===== RUN ANALYSIS =====
tryCatch({
  RnBeads::rnb.run.analysis(
    dir.reports = report.dir,
    sample.sheet = sample.annotation,
    data.dir = idat.dir,
    data.type = "idat.dir"
  )
}, error = function(e) {
  message("Analysis failed with error: ", e$message)
  stop("RnBeads pipeline aborted.")
})