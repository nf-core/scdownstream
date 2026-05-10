#!/usr/bin/env Rscript

options(future.globals.maxSize = 8 * 1024^3)

library(future)
library(Seurat)
library(anndataR)

set.seed(0)

plan("multicore", workers = max(1L, as.integer("${task.cpus}")))

adata <- read_h5ad("${h5ad}")

seurat_obj <- adata\$as_Seurat(x_mapping = "counts")

# Split layers by batch — V5 multi-layer integration workflow
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj@meta.data[["${batch_col}"]])

# SCTransform across all batches at once; override variable feature selection
# to use all genes since HVG selection was already performed upstream
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
VariableFeatures(seurat_obj) <- rownames(seurat_obj)

seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Correct the PCA latent space across batches (analogous to scVI latent correction)
seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    method = CCAIntegration,
    normalization.method = "SCT",
    new.reduction = "X_emb",
    verbose = FALSE
)

# Extract embeddings and align to adata cell order
emb <- Embeddings(seurat_obj, reduction = "X_emb")
emb <- emb[match(rownames(adata\$obs), rownames(emb)), ]
rownames(emb) <- NULL
colnames(emb) <- NULL
emb <- round(emb, 10)

adata\$obsm\$X_emb <- emb

write_h5ad(adata, "${prefix}.h5ad")

################################################
## VERSIONS FILE                              ##
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
