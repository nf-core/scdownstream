#!/usr/bin/env python3

import importlib.metadata
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
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

adata = read_anndata("${anndata}")
adata_proc = adata.copy()
prefix = "${prefix}"
batch_col = "${batch_col}"
input_layer = "${input_layer}"
log_normalize = "${log_normalize}" == "true"

if input_layer != "X" and input_layer not in adata_proc.layers:
    raise ValueError(
        f"input_layer {input_layer!r} is not present in adata.layers (available: {list(adata_proc.layers.keys())})"
    )

if input_layer != "X":
    adata_proc.X = adata_proc.layers[input_layer]

if log_normalize:
    target_sum = float(np.median(np.asarray(adata_proc.X.sum(axis=1)).ravel()))
    sc.pp.normalize_total(adata_proc, target_sum=target_sum)
    sc.pp.log1p(adata_proc)

sc.pp.scale(adata_proc, max_value=10)
sc.pp.pca(adata_proc, n_comps=30, zero_center=False)

sce.pp.scanorama_integrate(adata_proc, batch_col)

adata.obsm["X_scanorama"] = adata_proc.obsm["X_scanorama"]
adata.obsm["X_emb"] = adata_proc.obsm["X_scanorama"]
write_anndata(adata, prefix, "${anndata_format}")

pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names).to_pickle(f"X_{prefix}.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "scanorama": importlib.metadata.version("scanorama"),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
