#!/usr/bin/env Rscript

library(scDblFinder)
library(tidyverse)
library(SingleCellExperiment)
library(BiocParallel)
library(anndataR)

adata <- read_h5ad("${h5ad}")
sce <- adata\$as_SingleCellExperiment()

# Set the param to a specified RNG seed for reproducibility
bp <- MulticoreParam(workers = multicoreWorkers(), RNGseed=123)


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
multiplet_rate <- multiplet_rates_10x %>%
  dplyr::filter(Recovered_cells < ncol(sce)) %>%
  dplyr::slice(which.max(Recovered_cells)) %>%
  dplyr::pull(Multiplet_rate) %>%
  as.numeric()

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

# Restore original cell names
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

# Convert back to AnnData and save
adata_processed <- as_AnnData(sce)
write_h5ad(adata_processed, "${prefix}.h5ad")

# Extract predictions for doublet removal step
# Create a binary doublet call based on class
# Ensure we have valid row names
if (is.null(colnames(sce)) || length(colnames(sce)) != ncol(sce)) {
    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
}

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
