#!/usr/bin/env Rscript

library(scDblFinder)
library(tidyverse)
library(SingleCellExperiment)
library(BiocParallel)
library(anndataR)

# Set random seed for reproducibility
set.seed(123)

adata <- read_h5ad("${h5ad}")
sce <- adata\$as_SingleCellExperiment()

num_threads <- max(1L, as.integer("${task.cpus}"))
bp <- MulticoreParam(workers = num_threads, RNGseed = 123)

# Parse per-sample doublet rate from Nextflow input. If unavailable, let
# scDblFinder estimate dbr internally (recommended default for 10X data).
dbr_raw <- trimws("${dbr}")
dbr <- suppressWarnings(as.numeric(dbr_raw))

# Fetch sample batch information from nextflow metadata
batch_col <- trimws("${batch_col ?: ''}")

# Initialize sample groups(batches) variable to NULL to prevent object not found error
sample_groups <- NULL

# Check that the specified batch information exists in the analysis object
if (nzchar(batch_col)) {
  if (!(batch_col %in% colnames(colData(sce)))) {
    stop(
      "Batch column '",
      batch_col,
      "' was requested for scDblFinder samples but is not present in the input data. Available columns: ",
      paste(colnames(colData(sce)), collapse = ", ")
    )
  }

# Check that the batch column does not contain NAs
  sample_groups <- colData(sce)[[batch_col]]
  if (any(is.na(sample_groups))) {
    stop("Batch column '", batch_col, "' contains NA values; cannot split scDblFinder by sample.")
  }

# Assign the batch column from as the 'samples'
  sample_groups <- as.vector(sample_groups)
  message(paste0("Using batch column for scDblFinder samples: ", batch_col))
}
# Check for presence of metadata (doublet rate and batch annotation for sample groups)
if (is.na(dbr)) {
  message("No valid doublet_rate provided; using scDblFinder internal dbr estimation")
  dbr <- NULL
} else {
  message(paste0("Using provided doublet_rate (dbr): ", dbr))
}

scdblfinder_args <- list(
  assays(sce)[[1]],
  BPPARAM = bp,
  dbr = dbr
)

if (!is.null(sample_groups)) {
  scdblfinder_args\$samples <- sample_groups
}

# Run scDblFinder on the counts matrix (first assay)
# scDblFinder creates artificial doublets internally and returns a new SCE object
sce <- do.call(scDblFinder, scdblfinder_args)

# Generate a summary table
message("scDblFinder results summary:")
print(table(sce\$scDblFinder.class))

# Rename scDblFinder.* columns for consistency with other doublet methods.
# Replace prefix first, then replace any remaining dots with underscores.
idx <- grep("^scDblFinder\\\\.", colnames(colData(sce)))
colnames(colData(sce))[idx] <- gsub(
  "\\\\.",
  "_",
  sub("^scDblFinder\\\\.", "scdblfinder_", colnames(colData(sce))[idx])
)

# The doublet calls must stay keyed by the original cell barcodes. If they are not
# present here, something went wrong during conversion or scDblFinder processing and
# we should fail instead of inventing replacement identifiers.
if (is.null(colnames(sce)) || length(colnames(sce)) != ncol(sce)) {
  stop("scDblFinder output is missing valid cell barcodes; cannot write aligned h5ad and prediction outputs.")
}

# Write the updated SingleCellExperiment directly as h5ad.
write_h5ad(sce, "${prefix}.h5ad")

# Extract predictions for doublet removal step
# Create a binary doublet call based on class

# Create predictions vector
doublet_calls <- colData(sce)\$scdblfinder_class == "doublet"

# Create data frame without row.names first, then add them
predictions <- data.frame(doublet = doublet_calls)
row.names(predictions) <- colnames(sce)

colnames(predictions) <- "${prefix}"

# Save predictions to CSV
write.csv(predictions, "${prefix}.csv")

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
scDblFinder.version <- as.character(packageVersion('scDblFinder'))

writeLines(
    c(
        '"${task.process}":',
        paste('    R:', r.version),
        paste('    scDblFinder:', scDblFinder.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
