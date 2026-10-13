#!/usr/bin/env python3

import importlib.metadata
import platform

import anndata as ad
import yaml


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        return ad.read_zarr(path)
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(f"{stem}.zarr")


column = "$column"
adata = read_anndata("${anndata}")

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
write_anndata(adata_subset, "${prefix}", "${anndata_format}")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "numpy": importlib.metadata.version("numpy"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
