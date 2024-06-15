#!/usr/bin/env Rscript

library(Seurat)

sce <- readRDS("${rds}")

seurat_obj <- as.Seurat(sce, data = NULL)
seurat_obj <- seurat_obj[, unname(which(colSums(GetAssayData(seurat_obj)) != 0))]

batches <- SplitObject(seurat_obj, split.by = "batch")
batches <- lapply(X = batches, FUN = SCTransform, assay = "originalexp")

features <- SelectIntegrationFeatures(object.list = batches, nfeatures = 2000)
batches <- PrepSCTIntegration(object.list = batches, anchor.features = features)
batches <- lapply(X = batches, FUN = RunPCA, features = features, verbose=F)

anchors <- FindIntegrationAnchors(
    object.list = batches,
    normalization.method = "SCT",
    anchor.features = features
)

integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT"
)

saveRDS(integrated, "${prefix}.rds")

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
seurat.version <- as.character(packageVersion('Seurat'))

writeLines(
    c(
        '"${task.process}":',
        paste('    R:', r.version),
        paste('    Seurat:', seurat.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
