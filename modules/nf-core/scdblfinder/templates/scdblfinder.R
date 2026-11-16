#!/usr/bin/env Rscript

library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(anndataR)
library(arrow)

set.seed(123)

sce <- read_h5ad("${h5ad}")\$as_SingleCellExperiment()

# Extra arguments for scDblFinder from task.ext.args, written as R named
# arguments, e.g. dbr = 0.08, samples = "batch", clusters = TRUE
extra <- list(${args})

# A single string given as 'samples' is read as the name of a column in the
# cell metadata (obs), so the cells can be split by sample or capture
if (is.character(extra\$samples) && length(extra\$samples) == 1) {
    batch_col <- extra\$samples
    if (!(batch_col %in% colnames(colData(sce)))) {
        stop(
            "Column '", batch_col, "' is not in the cell metadata. Available columns: ",
            paste(colnames(colData(sce)), collapse = ", ")
        )
    }
    extra\$samples <- as.vector(colData(sce)[[batch_col]])
    if (anyNA(extra\$samples)) {
        stop("Column '", batch_col, "' contains NA values.")
    }
}

workers <- max(1L, as.integer("${task.cpus}"))
res <- do.call(
    scDblFinder,
    c(
        list(assay(sce, 1), BPPARAM = MulticoreParam(workers = workers, RNGseed = 123)),
        extra
    )
)

# Add the scDblFinder results to the original object so the other
# cell metadata is kept
cols <- grep("^scDblFinder", colnames(colData(res)), value = TRUE)
for (nm in cols) {
    sce[[nm]] <- res[[nm]]
}

message("scDblFinder results:")
print(table(sce\$scDblFinder.class))

write_h5ad(sce, "${prefix}.h5ad", compression = "gzip")

cell_ids <- colnames(sce)
if (is.null(cell_ids) || anyNA(cell_ids) || anyDuplicated(cell_ids) || !is.character(cell_ids)) {
    stop("Doublet predictions require unique, non-null string cell IDs")
}
scores <- as.data.frame(colData(sce))[, cols, drop = FALSE]
scores <- cbind(data.frame(`__index__` = cell_ids, check.names = FALSE), scores)
write_parquet(scores, "${prefix}.scores.parquet", compression = "snappy")

calls <- data.frame(sce\$scDblFinder.class == "doublet", row.names = colnames(sce))
colnames(calls) <- "${prefix}"
calls <- cbind(data.frame(`__index__` = cell_ids, check.names = FALSE), calls)
write_parquet(calls, "${prefix}.predictions.parquet", compression = "snappy")

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    scdblfinder:', as.character(packageVersion('scDblFinder'))),
        paste('    anndatar:', as.character(packageVersion('anndataR'))),
        paste('    singlecellexperiment:', as.character(packageVersion('SingleCellExperiment'))),
        paste('    biocparallel:', as.character(packageVersion('BiocParallel'))),
        paste('    rhdf5:', as.character(packageVersion('rhdf5'))),
        paste('    arrow:', as.character(packageVersion('arrow')))
    ),
    'versions.yml'
)
