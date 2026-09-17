#!/usr/bin/env Rscript

library(anndataR)

counts_layer <- "${counts_layer}"
prefix <- "${prefix}"
h5ad <- "${h5ad}"

adata <- read_h5ad(h5ad)
sce <- adata\$as_SingleCellExperiment()

SummarizedExperiment::assayNames(sce)[SummarizedExperiment::assayNames(sce) == counts_layer] <- "counts"

saveRDS(sce, paste0(prefix, ".rds"))

ver <- packageVersion("anndataR")

writeLines(
    c(
        "${task.process}:",
        paste0("    anndataR: ", ver)
    ),
    "versions.yml"
)
