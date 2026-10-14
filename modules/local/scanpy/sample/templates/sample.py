#!/usr/bin/env python3

import platform
from pathlib import Path

import anndata as ad
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits


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


threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


adata = read_anndata("${anndata}")
prefix = "${prefix}"
n = "${n ?: ''}"
fraction = "${fraction ?: ''}"

kwargs = {"rng": 0}

if n:
    kwargs["n"] = int(n)
elif fraction:
    kwargs["fraction"] = float(fraction)
else:
    raise ValueError("Either n or fraction must be set")

if "n" in kwargs and kwargs["n"] >= adata.n_obs:
    print(f"Warning: n is greater than the number of cells in the dataset ({adata.n_obs}). Using the entire dataset.")
else:
    sc.pp.sample(adata, **kwargs)

write_anndata(adata, prefix, "${anndata_format}")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
