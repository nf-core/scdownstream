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

top_idx <- ord[seq_len(n_keep)]

adata_out <- read_h5ad("${h5ad}")
if (!("counts" %in% names(adata_out\$layers))) {
    adata_out\$layers[["counts"]] <- adata_out\$X
}

var_idx <- match(rownames(sce), adata_out\$var_names)
if (any(is.na(var_idx))) {
    missing <- rownames(sce)[is.na(var_idx)]
    stop(
        "SCE genes missing from input AnnData var_names: ",
        paste(head(missing, 5), collapse = ", ")
    )
}
adata_out <- adata_out[, var_idx]

sce_genes <- rownames(sce)[top_idx]
adata_genes <- adata_out\$var_names[top_idx]
if (!identical(sce_genes, adata_genes)) {
    stop(
        "Gene order mismatch between deviance selection and input AnnData at selected indices. ",
        "First mismatch: ",
        paste(head(sce_genes[sce_genes != adata_genes], 5), collapse = ", ")
    )
}

obs_meta <- adata_out\$obs
adata_out <- adata_out[, top_idx]
adata_out <- adata_out\$as_InMemoryAnnData()
# HDF5-backed views can drop obs when row names are duplicated; restore metadata when possible.
if (ncol(adata_out\$obs) == 0 && ncol(obs_meta) > 0) {
    adata_out\$obs <- obs_meta
}
adata_out\$var[["binomial_deviance"]] <- deviance[top_idx]
adata_out\$var[["highly_deviant"]] <- highly_deviant[top_idx]
write_h5ad(adata_out, "${prefix}.h5ad")

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
scry.version <- as.character(packageVersion('scry'))

writeLines(
    c(
        '"${task.process}":',
        paste('    R:', r.version),
        paste('    scry:', scry.version)
    ),
'versions.yml')
