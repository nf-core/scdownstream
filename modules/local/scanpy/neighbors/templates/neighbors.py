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


adata = read_anndata("${anndata}", backed="r")
prefix = "${prefix}"

kwargs = {"use_rep": "${rep}"}
n_pcs = "${n_pcs ?: ''}"
if n_pcs:
    kwargs["n_pcs"] = int(n_pcs)

sc.pp.neighbors(adata, **kwargs)

write_anndata(adata, prefix, "${anndata_format}")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
