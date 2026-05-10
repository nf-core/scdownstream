#!/usr/bin/env Rscript

library(SoupX)
library(anndataR)
library(Seurat)

# Read the AnnData objects and convert directly to Seurat
adata <- read_h5ad("${h5ad}")
seu <- adata\$as_Seurat(x_mapping = "counts")

# Read the raw AnnData object
adata_raw <- read_h5ad("${raw}")
seu_raw <- adata_raw\$as_Seurat(x_mapping = "counts")

# Get layer to use (default to "X" if not specified)
use_layer <- "${input_layer}"
if (use_layer != "X") {
    if (use_layer %in% names(adata\$layers)) {
        # Replace counts in the Seurat object with the specified layer
        counts_matrix <- t(adata\$layers[[use_layer]])
        rownames(counts_matrix) <- rownames(seu)
        colnames(counts_matrix) <- colnames(seu)
        seu[["RNA"]] <- CreateAssayObject(counts = counts_matrix)
    } else {
        stop(paste("Layer", use_layer, "not found in AnnData object"))
    }
}

# Preprocessing with Seurat workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu)
seu <- FindClusters(seu, resolution = ${cluster_resolution})
soupx_groups <- Idents(seu)

# Create SoupChannel
data <- GetAssayData(seu, layer = "counts")
data_raw <- GetAssayData(seu_raw, layer = "counts")
sc <- SoupChannel(data_raw, data, calcSoupProfile = FALSE)

# Set soup profile - following SoupX vignette approach
soupProf <- data.frame(
    row.names = rownames(data),
    est = rowSums(data)/sum(data),
    counts = rowSums(data)
)
sc <- setSoupProfile(sc, soupProf)

# Set clusters and estimate contamination
sc <- setClusters(sc, soupx_groups)
sc <- autoEstCont(sc, doPlot = FALSE)

# Adjust counts
out <- adjustCounts(sc, roundToInt = FALSE)

# Update the original AnnData object with ambient layer
if ("${output_layer}" == "X") {
    adata\$X <- t(out)
} else {
    adata\$layers[["${output_layer}"]] <- t(out)
}

# Save the output
write_h5ad(adata, "${prefix}.h5ad")

# Write version information
writeLines(
    c(
        '"${task.process}":',
        paste('    r:', paste(R.Version()\$major, R.Version()\$minor, sep = ".")),
        paste('    anndataR:', as.character(packageVersion('anndataR'))),
        paste('    SoupX:', as.character(packageVersion('SoupX'))),
        paste('    Seurat:', as.character(packageVersion('Seurat')))
    ),
    'versions.yml'
)
