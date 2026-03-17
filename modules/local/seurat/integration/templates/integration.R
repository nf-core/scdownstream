#!/usr/bin/env Rscript

options(future.globals.maxSize = 8 * 1024^3)

library(Seurat)
library(anndataR)

# Set random seed
set.seed(0)

adata <- read_h5ad("${h5ad}")

seurat_obj <- adata\$as_Seurat()

batches <- SplitObject(seurat_obj, split.by = "${batch_col}")
batches <- lapply(X = batches, FUN = SCTransform, assay = "RNA")

features <- SelectIntegrationFeatures(object.list = batches, nfeatures = nrow(seurat_obj))
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

integrated <- RunPCA(integrated, reduction.name = "X_emb")

# Extract embeddings
emb <- Embeddings(integrated, reduction = "X_emb")

# Make sure the emb row order matches
emb <- emb[match(rownames(adata\$obs), rownames(emb)), ]

# Remove row and column names from the matrix
rownames(emb) <- NULL
colnames(emb) <- NULL

# Round to 10 decimal places
emb <- round(emb, 10)

# Add embeddings to adata
adata\$obsm\$X_emb <- emb

write_h5ad(adata, "${prefix}.h5ad")

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
seurat.version <- as.character(packageVersion('Seurat'))
anndataR.version <- as.character(packageVersion('anndataR'))

writeLines(
    c(
        '"${task.process}":',
        paste('    R:', r.version),
        paste('    Seurat:', seurat.version),
        paste('    anndataR:', anndataR.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
