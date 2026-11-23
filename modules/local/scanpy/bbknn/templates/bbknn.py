#!/usr/bin/env python3

import platform
from importlib.metadata import version

import scanpy as sc
import yaml
from anndata.acc import A
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
input_layer = "${input_layer}"
log_normalize = "${log_normalize}" == "true"

if input_layer != "X" and input_layer not in adata.layers:
    raise ValueError(
        f"input_layer {input_layer!r} is not present in adata.layers (available: {list(adata.layers.keys())})"
    )

analysis_layer = None if input_layer == "X" else input_layer
if log_normalize:
    output_layer = "bbknn_log1p"
    if output_layer in adata.layers:
        raise ValueError(f"Normalisation output layer {output_layer!r} already exists")
    source = adata.X if analysis_layer is None else adata.layers[analysis_layer]
    adata.layers[output_layer] = source.copy()
    sc.pp.normalize_total(adata, layer=output_layer)
    sc.pp.log1p(adata, layer=output_layer)
    analysis_layer = output_layer

sc.tl.pca(adata, layer=analysis_layer)

sc.pp.bbknn(
    adata,
    neighbors_within_batch=25 if adata.n_obs >= 100_000 else 3,
    batch=A.obs["${batch_col}"],
    use=A.obsm["X_pca"],
    rng=0,
)

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

versions = {
    "${task.process}": {
        "anndata": version("anndata"),
        "pandas": version("pandas"),
        "python": platform.python_version(),
        "scanpy": version("scanpy"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
