#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# bioconda installs celldex in a post-link step that pixi-based Wave images skip,
# so its fetchReference() backends gypsum and alabaster are used directly.
library(gypsum)
library(alabaster.base)
library(alabaster.se)
library(SingleCellExperiment)
library(yaml)
library(HDF5Array)

prefix <- "${prefix}"
ref_name <- "${ref}"
ref_version <- "${version}"

print(paste("Attempting to fetch reference:", ref_name, ref_version))

# Equivalent to celldex::fetchReference(ref_name, ref_version, cache = "./")
version_path <- saveVersion("celldex", ref_name, ref_version, cache = "./")
reference <- readObject(version_path)

# Save SummarizedExperiment to HDF5 files
saveHDF5SummarizedExperiment(
  reference,
  dir = prefix,
  replace = TRUE
)
# Compress the HDF5 files into a tar archive
tar(tarfile = paste0(prefix, ".tar"), files = prefix)

versions <- list(
  "${task.process}" = list(
    R = R.version.string,
    gypsum = as.character(packageVersion("gypsum")),
    alabaster.base = as.character(packageVersion("alabaster.base")),
    alabaster.se = as.character(packageVersion("alabaster.se")),
    yaml = as.character(packageVersion("yaml")),
    SingleCellExperiment = as.character(
      packageVersion("SingleCellExperiment")
    ),
    HDF5Array = as.character(packageVersion("HDF5Array"))
  )
)
write_yaml(x = versions, file = "versions.yml")
