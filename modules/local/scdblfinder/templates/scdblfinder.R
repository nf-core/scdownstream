#!/usr/bin/env Rscript

library(scDblFinder)
library(tidyverse)
library(SingleCellExperiment)
library(BiocParallel)
library(anndataR)

adata <- read_h5ad("${h5ad}")
sce <- adata$as_SingleCellExperiment()

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

# Run scDblFinder on the REAL data (not mock data!)
# scDblFinder creates artificial doublets internally
set.seed(123)
sce <- scDblFinder(
    sce,
    BPPARAM = bp,
    dbr = multiplet_rate,
    artificialDoublets = ncol(sce)
)

# Generate a summary table
message("scDblFinder results summary:")
print(table(sce\$scDblFinder.class))

# Rename scDblFinder.* columns for consistency with other doublet methods
scdbl_cols <- grep("^scDblFinder\\\\.", colnames(colData(sce)), value = TRUE)
new_scdbl_cols <- paste0("scdblfinder_", gsub("^scDblFinder\\\\.", "", gsub("\\\\.", "_", scdbl_cols)))

# Rename columns in colData(sce)
for (i in seq_along(scdbl_cols)) {
  colData(sce)[[new_scdbl_cols[i]]] <- colData(sce)[[scdbl_cols[i]]]
  colData(sce)[[scdbl_cols[i]]] <- NULL  # Remove the original column
}

# Convert back to AnnData and save
adata_processed <- as_AnnData(sce)
write_h5ad(adata_processed, "${prefix}.h5ad")

# Extract predictions for doublet removal step
# Create a binary doublet call based on class
predictions <- data.frame(
    doublet = colData(sce)\$scdblfinder_class == "doublet",
    row.names = colnames(sce)
)
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
