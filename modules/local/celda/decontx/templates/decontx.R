#!/usr/bin/env Rscript

library(anndataR)
library(celda)

# Read the AnnData object
adata <- read_h5ad("${h5ad}")

# Convert to SingleCellExperiment, mapping only the target input to "counts"
# assays_mapping = FALSE avoids duplicate name conflicts when a "counts" layer already exists
input_layer <- "${input_layer}"
if (input_layer == "X") {
    sce <- adata\$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)
} else {
    sce <- adata\$as_SingleCellExperiment(x_mapping = FALSE, assays_mapping = c(counts = input_layer))
}

# Prepare parameters for decontX
params <- list()
params\$assayName <- "counts"

# Handle batch information if available
batch_col <- "${batch_col ?: ''}"

if (batch_col != "") {
    # Plausibility check 1: Does the column exist?
    if (!(batch_col %in% colnames(adata\$obs))) {
        cat("ERROR: Batch column '", batch_col, "' not found in obs columns.\\n", sep="")
        cat("Available columns in obs:", paste(colnames(adata\$obs), collapse=", "), "\\n")
        stop("Batch column '", batch_col, "' does not exist in the AnnData object. Available columns: ", paste(colnames(adata\$obs), collapse=", "))
    }

    batch_values <- adata\$obs[[batch_col]]

    # Plausibility check 2: Does it contain NA/NaN values?
    na_count <- sum(is.na(batch_values))
    if (na_count > 0) {
        unique_vals <- unique(batch_values)
        cat("ERROR: Batch column '", batch_col, "' contains ", na_count, " NA/NaN value(s).\\n", sep="")
        cat("Unique batch values (including NA/NaN):", paste(unique_vals, collapse=", "), "\\n")
        stop("Batch column '", batch_col, "' contains NA/NaN values. All batch values must be non-NA. Unique values: ", paste(unique_vals, collapse=", "))
    }

    # Both checks passed - use the batch column
    params\$batch <- batch_values
}

# Handle background data if available
raw_path <- "${raw}"
if (file.exists(raw_path)) {
    raw <- read_h5ad(raw_path)
    raw_sce <- raw\$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)
    params\$background <- raw_sce

    # If a batch column is provided, we also need to add it to the background data
    if (!is.null(params\$batch)) {
        params\$bgBatch <- raw\$obs[[batch_col]]
    }
}

# Run decontX with parameters
corrected <- do.call(decontX, c(list(sce), params))

adata_corrected <- as_AnnData(corrected)

# Convert back to AnnData and update layers
if ("${output_layer}" == "X") {
    adata\$X <- adata_corrected\$layers[["decontXcounts"]]
} else {
    adata\$layers[["${output_layer}"]] <- adata_corrected\$layers[["decontXcounts"]]
}

# Save the output
write_h5ad(adata, "${prefix}.h5ad")

# Write version information
writeLines(
    c(
        '"${task.process}":',
        paste('    r:', paste(version\$major, version\$minor, sep = ".")),
        paste('    anndataR:', as.character(packageVersion('anndataR'))),
        paste('    celda:', as.character(packageVersion('celda')))
    ),
'versions.yml')
