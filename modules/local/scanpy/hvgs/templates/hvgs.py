#!/usr/bin/env python3

import platform

import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
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

if n_hvgs > 0 and adata.n_vars <= n_hvgs:
    # When available genes <= positive requested limit, select all genes
    adata.var["highly_variable"] = True
    adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")
elif adata.n_vars > n_hvgs:
    kwargs = {}

    if batch_key:
        kwargs["batch_key"] = batch_key

    # If an actual limit is provided, use it
    # Otherwise, scanpy will automatically determine the number of highly variable genes
    if n_hvgs > 0:
        kwargs["n_top_genes"] = n_hvgs

    if input_layer != "X" and input_layer not in adata.layers:
        raise ValueError(
            f"input_layer {input_layer!r} is not present in adata.layers (available: {list(adata.layers.keys())})"
        )

    analysis_layer = None if input_layer == "X" else input_layer
    if log_normalize:
        output_layer = "hvgs_log1p"
        if output_layer in adata.layers:
            raise ValueError(f"Normalisation output layer {output_layer!r} already exists")
        source = adata.X if analysis_layer is None else adata.layers[analysis_layer]
        adata.layers[output_layer] = source.copy()
        sc.pp.normalize_total(adata, layer=output_layer, target_sum=None)
        sc.pp.log1p(adata, layer=output_layer)
        analysis_layer = output_layer

    sc.pp.highly_variable_genes(adata, layer=analysis_layer, **kwargs)

    adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")

    adata = adata[:, adata.var["highly_variable"]]

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
