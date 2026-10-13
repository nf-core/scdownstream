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


adata = read_anndata("${anndata}")

axis = "$axis"
column = "$column"

df = adata.obs if axis == "obs" else adata.var

assert column in df.columns, f"Column '{column}' not found in {axis}"
df.index = df[column]
df.index.name = None

if axis == "obs":
    adata.obs = df
else:
    adata.var = df

write_anndata(adata, "${prefix}", "${anndata_format}")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
