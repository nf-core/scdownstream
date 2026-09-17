#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

library(celldex)
library(SingleCellExperiment)
library(yaml)
library(HDF5Array)

# Split the reference into ref_name and ref_version based on __
prefix <- "${prefix}"
ref_name <- "${ref}"
ref_version <- "${version}"

print(paste("Attempting to fetch reference:", ref_name, ref_version))

reference <- fetchReference(ref_name, ref_version, cache = "./")

# Save SummarizedExperiment to HDF5 files
saveHDF5SummarizedExperiment(
  reference,
  dir = prefix,
  replace = TRUE
)
# Compress the HDF5 files into a tar archive
tar(tarfile = paste0(prefix, ".tar"), files = prefix)

# Capturing version information, as before
versions <- list(
  "${task.process}" = list(
    R = R.version.string,
    celldex = as.character(packageVersion("celldex")),
    yaml = as.character(packageVersion("yaml")),
    SingleCellExperiment = as.character(
      packageVersion("SingleCellExperiment")
    ),
    HDF5Array = as.character(packageVersion("HDF5Array"))
  )
)
# Write versions info into a YAML file, as before
write_yaml(x = versions, file = "versions.yml")
