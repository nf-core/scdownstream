#!/usr/bin/env python3

import os

os.environ["KMP_AFFINITY"] = "disabled"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform

import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
n_genes = int("${n_genes}")
batch_key = "${batch_key}"
counts_layer = "${counts_layer}"

if n_genes < 0:
    raise ValueError(
        "pearson_residuals_hvgs requires integration_n_features >= 0; "
        "negative values are not supported for this feature selection method."
    )

if n_genes <= 0:
    n_genes = 4000

if "${excluded_genes}":
    with open("${excluded_genes}") as f:
        excluded_genes = [line.strip() for line in f if line.strip()]
    mask = ~adata.var_names.isin(excluded_genes)
    adata = adata[:, mask].copy()

if counts_layer == "X":
    raw_counts = adata.X.copy()
    pearson_layer = None
else:
    if counts_layer not in adata.layers:
        raise ValueError(
            f"counts_layer '{counts_layer}' was requested but is not present in adata.layers"
        )
    raw_counts = adata.layers[counts_layer].copy()
    pearson_layer = counts_layer

if adata.n_vars > n_genes:
    kwargs = {
        "flavor": "pearson_residuals",
        "n_top_genes": n_genes,
        "chunksize": 1000,
    }

    if pearson_layer:
        kwargs["layer"] = pearson_layer

    if batch_key:
        kwargs["batch_key"] = batch_key

    sc.experimental.pp.highly_variable_genes(adata, **kwargs)

    adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")

    adata.X = raw_counts
    adata = adata[:, adata.var["highly_variable"]]
else:
    adata.var["highly_variable"] = True
    adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")
    adata.X = raw_counts

adata.write_h5ad(f"{prefix}.h5ad")

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
