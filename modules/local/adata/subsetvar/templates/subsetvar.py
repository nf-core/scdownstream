#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform

import anndata as ad
import yaml

column = "$column"
adata = ad.read_h5ad("$h5ad")

if column not in adata.var.columns:
    raise ValueError(f"Column '{column}' not found in adata.var")

mask = adata.var[column]
if mask.dtype != bool:
    raise ValueError(f"Column '{column}' must have boolean dtype, but got {mask.dtype}")

if mask.isna().any():
    raise ValueError(f"Column '{column}' contains null values")

selected = mask.to_numpy()
if not selected.any():
    raise ValueError(f"Column '{column}' selects zero genes")

adata_subset = adata[:, selected].copy()
ad.settings.allow_write_nullable_strings = True
adata_subset.write_h5ad("${prefix}.h5ad")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "numpy": importlib.metadata.version("numpy"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
