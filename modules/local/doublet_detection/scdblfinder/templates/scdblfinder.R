#!/usr/bin/env Rscript

library(scDblFinder)
library(tidyverse)
library(SingleCellExperiment)
library(BiocParallel)
library(anndataR)

adata <- read_h5ad("${h5ad}")
sce <- adata\$as_SingleCellExperiment()

# Set the param to a specified RNG seed for reproducibility
num_threads <- max(1L, as.integer("${task.cpus}"))
bp <- MulticoreParam(workers = num_threads, RNGseed = 123)

# 10 Genomics Doublet Rate calculator used to get multiplet rate if not provided
# 10X multiplet rate table(https://rpubs.com/kenneditodd/doublet_finder_example)
multiplet_rates_10x <- data.frame(
  "Multiplet_rate" = c(0.004, 0.008, 0.0160, 0.023, 0.031,
                        0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
  "Loaded_cells" = c(800, 1600, 3200, 4800, 6400, 8000, 9600,
                     11200, 12800, 14400, 16000),
  "Recovered_cells" = c(500, 1000, 2000, 3000, 4000, 5000, 6000,
                        7000, 8000, 9000, 10000)
)

# Adjust to use the number of cells in the SCE object
idx <- findInterval(ncol(sce), multiplet_rates_10x\$Recovered_cells)
if (idx < 1L) idx <- 1L
if (idx > nrow(multiplet_rates_10x)) idx <- nrow(multiplet_rates_10x)

multiplet_rate <- as.numeric(multiplet_rates_10x\$Multiplet_rate[idx])
message(paste0("Setting multiplet rate to ", multiplet_rate, " for ", ncol(sce), " cells"))

# Save original cell names and count before overwriting sce
original_cell_names <- colnames(sce)
n_cells <- ncol(sce)

# Run scDblFinder on the counts matrix (first assay)
# scDblFinder creates artificial doublets internally and returns a new SCE
set.seed(123)
sce <- scDblFinder(
    assays(sce)[[1]],
    BPPARAM = bp,
    dbr = multiplet_rate,
    artificialDoublets = n_cells
)

# Restore the input barcodes because running scDblFinder on the just the assay matrix above can
# return a new SCE whose column names no longer match the original AnnData cell IDs.
# Keeping the original names is required so the output h5ad obs_names and CSV rows
# still map back to the same cells seen by downstream steps.
if (!is.null(original_cell_names) && length(original_cell_names) == ncol(sce)) {
    colnames(sce) <- original_cell_names
}

# Generate a summary table
message("scDblFinder results summary:")
print(table(sce\$scDblFinder.class))

# Rename scDblFinder.* columns for consistency with other doublet methods
scdbl_cols <- grep("^scDblFinder\\\\.", colnames(colData(sce)), value = TRUE)

# First remove "scDblFinder." prefix, THEN replace remaining dots with underscores
new_scdbl_cols <- paste0("scdblfinder_", gsub("\\\\.", "_", gsub("^scDblFinder\\\\.", "", scdbl_cols)))

# Rename columns in colData(sce) - create new columns first, then delete old ones
for (i in seq_along(scdbl_cols)) {
  colData(sce)[[new_scdbl_cols[i]]] <- colData(sce)[[scdbl_cols[i]]]
}
# Now delete old columns
for (col in scdbl_cols) {
  colData(sce)[[col]] <- NULL
}

# The doublet calls must stay keyed by the original cell barcodes. If they are not
# present here, something went wrong during conversion or scDblFinder processing and
# we should fail instead of inventing replacement identifiers.
if (is.null(colnames(sce)) || length(colnames(sce)) != ncol(sce)) {
  stop("scDblFinder output is missing valid cell barcodes; cannot write aligned h5ad and prediction outputs.")
}

# Write the updated SingleCellExperiment directly as h5ad, explicitly mapping the
# primary assay to AnnData X so downstream readers see a valid matrix field.
primary_assay <- assayNames(sce)[1]
if (is.na(primary_assay) || primary_assay == "") {
  stop("scDblFinder output is missing a primary assay; cannot write h5ad output.")
}
write_h5ad(sce, "${prefix}.h5ad", x_mapping = primary_assay)

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
