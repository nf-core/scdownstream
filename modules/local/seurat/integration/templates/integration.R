#!/usr/bin/env Rscript

library(Seurat)

sce <- readRDS("${rds}")

seurat_obj <- as.Seurat(sce)

batch_list <- SplitObject(seurat_obj, split.by = "batch")

anchors <- FindIntegrationAnchors(
  object.list = batch_list,
  anchor.features = 2000,
  scale = T,
  l2.norm = T,
  dims = 1:30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  eps = 0
)

integrated <- IntegrateData(
  anchorset = anchors,
  new.assay.name = "X_emb",
  features = NULL,
  features.to.integrate = NULL,
  dims = 1:30,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = F,
  do.cpp = T,
  eps = 0,
  verbose = T
)

saveRDS(integrated, "${prefix}.rds")