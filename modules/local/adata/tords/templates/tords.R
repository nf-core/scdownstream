#!/usr/bin/env Rscript

library(anndataR)

counts_layer <- "${counts_layer}"
prefix <- "${prefix}"
h5ad <- "${h5ad}"

adata <- read_h5ad(h5ad)
layer_names <- adata\$layers_keys()
if (counts_layer != "X" && !counts_layer %in% layer_names) {
    stop(sprintf(
        "Requested counts source '%s' is missing from '%s'. Select X or an available layer: %s",
        counts_layer, h5ad, paste(layer_names, collapse = ", ")
    ))
}
if (counts_layer == "X" && is.null(adata\$X)) {
    stop(sprintf("Requested counts source 'X' is missing from '%s'. Select an available counts layer.", h5ad))
}

assay_names <- make.unique(c("counts", layer_names))[-1]
assays_mapping <- setNames(layer_names, assay_names)
if (counts_layer == "X") {
    x_mapping <- "counts"
} else {
    assays_mapping <- c(counts = counts_layer, assays_mapping)
    x_mapping <- tail(make.unique(c(names(assays_mapping), "X")), 1)
}

# anndataR 1.0.2 treats x_mapping = FALSE as an assay name, not an exclusion.
sce <- adata\$as_SingleCellExperiment(
    x_mapping = x_mapping,
    assays_mapping = if (length(assays_mapping)) assays_mapping else FALSE
)
assay_names <- SummarizedExperiment::assayNames(sce)
SummarizedExperiment::assays(sce) <- SummarizedExperiment::assays(sce)[c("counts", setdiff(assay_names, "counts"))]

saveRDS(sce, paste0(prefix, ".rds"))

ver <- packageVersion("anndataR")

writeLines(
    c(
        "${task.process}:",
        paste0("    anndataR: ", ver)
    ),
    "versions.yml"
)
