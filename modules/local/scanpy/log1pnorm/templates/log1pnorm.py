#!/usr/bin/env python3

import platform
from pathlib import Path

import anndata as ad
import scanpy as sc
import yaml


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
prefix = "${prefix}"

if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

adata_norm = adata.copy()
adata_norm.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata_norm, target_sum=None)
sc.pp.log1p(adata_norm)

adata.layers["log1p"] = adata_norm.X.copy()

write_anndata(adata, prefix, "${anndata_format}")

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
