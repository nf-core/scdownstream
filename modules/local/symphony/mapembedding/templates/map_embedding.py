#!/usr/bin/env python3

import importlib.metadata
import platform
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import symphonypy as sp
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
adata_ref = sc.read_h5ad("reference/reference.h5ad")
prefix = "${prefix}"
batch_col = "${batch_col}"
counts_layer = "${counts_layer}"

if counts_layer != "X":
    adata_proc.X = adata_proc.layers[counts_layer]

target_sum = float(adata_ref.uns["normalize"]["target_sum"])
sc.pp.normalize_total(adata_proc, target_sum=target_sum)
sc.pp.log1p(adata_proc)

sp.tl.map_embedding(
    adata_proc,
    adata_ref,
    key=batch_col,
    transferred_adjusted_basis="X_symphony",
    use_genes_column="highly_variable",
)

adata.obsm["X_symphony"] = adata_proc.obsm["X_symphony"]
adata.obsm["X_emb"] = adata_proc.obsm["X_symphony"]

write_anndata(adata, prefix, "${anndata_format}")
pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names).to_pickle(f"X_{prefix}.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "symphonypy": importlib.metadata.version("symphonypy"),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
