#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
from threadpoolctl import threadpool_limits

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import yaml

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
n_hvgs = int("${n_hvgs}")
batch_key = "${batch_key}"

# Remove excluded genes from the anndata prior to identifying highly variable genes
if "${excluded_genes}":
    with open("${excluded_genes}", "r") as f:
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

    raw_counts = adata.X.copy()

    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, **kwargs)

    adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")

    adata.X = raw_counts
    adata = adata[:, adata.var["highly_variable"]]

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
