#!/usr/bin/env Rscript

library(scry)
library(SingleCellExperiment)
library(anndataR)

set.seed(123)

adata <- read_h5ad("${h5ad}")
sce <- adata\$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)

n_genes <- as.integer("${n_genes}")
if (n_genes <= 0) {
    n_genes <- 4000L
}

batch_col <- trimws("${batch_col}")

excluded_genes_path <- "${excluded_genes}"
if (nzchar(excluded_genes_path) && file.exists(excluded_genes_path)) {
    excluded_genes <- readLines(excluded_genes_path)
    excluded_genes <- excluded_genes[nzchar(trimws(excluded_genes))]
    if (length(excluded_genes) > 0) {
        keep <- !(rownames(sce) %in% excluded_genes)
        sce <- sce[keep, ]
    }
}

batch <- NULL
if (nzchar(batch_col)) {
    if (!(batch_col %in% colnames(colData(sce)))) {
        stop(
            "Batch column '",
            batch_col,
            "' was requested for deviance feature selection but is not present in the input data. Available columns: ",
            paste(colnames(colData(sce)), collapse = ", ")
        )
    }
    batch <- colData(sce)[[batch_col]]
    if (any(is.na(batch))) {
        stop("Batch column '", batch_col, "' contains NA values; cannot run batch-aware deviance feature selection.")
    }
}

sce <- devianceFeatureSelection(sce, assay = "counts", batch = batch)

deviance <- rowData(sce)\$binomial_deviance
ord <- order(deviance, decreasing = TRUE)
n_keep <- min(n_genes, nrow(sce))
top_genes <- rownames(sce)[ord[seq_len(n_keep)]]
highly_deviant <- rownames(sce) %in% top_genes

rowData(sce)\$binomial_deviance <- deviance
rowData(sce)\$highly_deviant <- highly_deviant

features_df <- data.frame(
    binomial_deviance = deviance,
    highly_deviant = highly_deviant,
    row.names = rownames(sce)
)
write.csv(features_df, "${prefix}_features.csv")

sce <- sce[top_genes, ]
write_h5ad(sce, "${prefix}.h5ad")

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
scry.version <- as.character(packageVersion('scry'))

writeLines(
    c(
        '"${task.process}":',
        paste('    R:', r.version),
        paste('    scry:', scry.version)
    ),
'versions.yml')
