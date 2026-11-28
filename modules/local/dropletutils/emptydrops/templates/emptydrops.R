#!/usr/bin/env Rscript

library(anndataR)
library(BiocParallel)
library(DropletUtils)
library(SingleCellExperiment)

adata <- read_h5ad("${h5ad}")
sce <- adata\$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)

lower <- as.numeric("${lower}")
fdr <- as.numeric("${fdr}")

set.seed(0)
result <- emptyDrops(
    counts(sce),
    lower = lower,
    BPPARAM = MulticoreParam(workers = ${task.cpus})
)

is_cell <- !is.na(result\$FDR) & result\$FDR <= fdr
cell_barcodes <- rownames(result)[is_cell]

if (length(cell_barcodes) == 0) {
    stop(
        "emptyDrops did not retain any barcodes at FDR <= ", fdr,
        " with lower = ", lower, ". Check the input matrix or adjust the thresholds."
    )
}

write.table(
    data.frame(barcode = cell_barcodes),
    "${prefix}_barcodes.csv",
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

write.csv(
    data.frame(barcode = rownames(result), as.data.frame(result)),
    "${prefix}_results.csv",
    row.names = FALSE
)

writeLines(
    c(
        '"${task.process}":',
        paste('    r:', paste(version\$major, version\$minor, sep = ".")),
        paste('    anndataR:', as.character(packageVersion('anndataR'))),
        paste('    DropletUtils:', as.character(packageVersion('DropletUtils')))
    ),
'versions.yml')
