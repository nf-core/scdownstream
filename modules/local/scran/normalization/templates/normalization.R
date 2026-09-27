#!/usr/bin/env Rscript

library(scran)
library(scater)
library(scuttle)
library(anndataR)
library(SingleCellExperiment)

adata <- read_h5ad("${h5ad}")

if (!("counts" %in% names(adata\$layers))) {
    adata\$layers[["counts"]] <- adata\$X
}

sce <- adata\$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)

lib_sizes <- librarySizeFactors(sce)
keep_cells <- lib_sizes > 0
if (!all(keep_cells)) {
    sce <- sce[, keep_cells]
}

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
size_factors <- sizeFactors(sce)
# Zero-library cells are removed above. Remaining invalid size factors are
# replaced with the smallest positive factor so every emitted cell keeps a
# valid normalisation scale and aligned output axes.
invalid <- !is.finite(size_factors) | size_factors <= 0
if (any(invalid)) {
    positive <- size_factors[!invalid]
    replacement <- if (length(positive) > 0) min(positive) else 1
    size_factors[invalid] <- replacement
    sizeFactors(sce) <- size_factors
}
sce <- logNormCounts(sce)

logcounts <- t(assay(sce, "logcounts"))
adata_out <- read_h5ad("${h5ad}")

if (!("counts" %in% names(adata_out\$layers))) {
    adata_out\$layers[["counts"]] <- adata_out\$X
}

if (!all(keep_cells)) {
    adata_out <- adata_out[keep_cells, ]
}

adata_out\$layers[["scran"]] <- logcounts
write_h5ad(adata_out, "${prefix}.h5ad")

r.version <- strsplit(version[["version.string"]], " ")[[1]][3]
scran.version <- as.character(packageVersion("scran"))

writeLines(
    c(
        '"${task.process}":',
        paste("    R:", r.version),
        paste("    scran:", scran.version)
    ),
    "versions.yml"
)
