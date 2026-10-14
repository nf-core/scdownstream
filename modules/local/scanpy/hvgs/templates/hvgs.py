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
n_hvgs = int("${n_hvgs}")
batch_key = "${batch_key}"
input_layer = "${input_layer}"
log_normalize = "${log_normalize}" == "true"

# Remove excluded genes from the anndata prior to identifying highly variable genes
if "${excluded_genes}":
    with open("${excluded_genes}") as f:
        excluded_genes = [line.strip() for line in f if line.strip()]
    mask = ~adata.var_names.isin(excluded_genes)
    adata = adata[:, mask].copy()

if adata.n_vars > n_hvgs:
    kwargs = {}

    if batch_key:
        kwargs["batch_key"] = batch_key

    # If an actual limit is provided, use it
    # Otherwise, scanpy will automatically determine the number of highly variable genes
    if n_hvgs > 0:
        kwargs["n_top_genes"] = n_hvgs

    raw_counts = adata.layers["counts"].copy() if "counts" in adata.layers else adata.X.copy()

    if input_layer != "X" and input_layer not in adata.layers:
        raise ValueError(
            f"input_layer {input_layer!r} is not present in adata.layers (available: {list(adata.layers.keys())})"
        )

    if input_layer != "X":
        adata.X = adata.layers[input_layer]

    if log_normalize:
        sc.pp.normalize_total(adata, target_sum=None)
        sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, **kwargs)

    adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")

    adata.X = raw_counts
    adata = adata[:, adata.var["highly_variable"]]

write_anndata(adata, prefix, "${anndata_format}")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
