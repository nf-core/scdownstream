#!/usr/bin/env python3

import platform
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy as sp
import yaml
from scipy.sparse import csc_matrix


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        # An absolute path keeps zarr from parsing colons in Nextflow prefixes as URL schemes
        return ad.read_zarr(Path(path).absolute())
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(Path(f"{stem}.zarr").absolute())


adata = read_anndata("${anndata}")

integration_methods = ["scanvi", "scvi", "symphony", "scimilarity", "seurat", "bbknn", "combat", "pca", "expimap"]
dimred_methods = ["umap", "tsne"]

for key in list(adata.obsm.keys()):
    if not any(key.endswith(dimred_method) for dimred_method in dimred_methods):
        del adata.obsm[key]

found_index = 1
for integration_method in integration_methods:
    for dimred_method in dimred_methods:
        key = f"X_{integration_method}-global_{dimred_method}"
        if key in adata.obsm.keys():
            new_key = f"X_{found_index}_{integration_method}-global_{dimred_method}"
            adata.obsm[new_key] = adata.obsm.pop(key)
            found_index += 1

# Delete everything in uns and layers (CELLxGENE can't display these anyway)
adata.uns = {}
adata.layers = {}

# Convert all float64 columns to float32
for df in [adata.obs, adata.var]:
    for col in df.columns:
        if df[col].dtype == np.float64:
            df[col] = df[col].astype(np.float32)

adata.X = csc_matrix(adata.X).astype(np.float32)
sc.pp.log1p(adata)

write_anndata(adata, "${prefix}", "${anndata_format}")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scipy": sp.__version__,
        "numpy": np.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
