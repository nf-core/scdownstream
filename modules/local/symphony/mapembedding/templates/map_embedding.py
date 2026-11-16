#!/usr/bin/env python3

import importlib.metadata
import platform
from pathlib import Path

import pandas as pd
import pyarrow as pa
import scanpy as sc
import symphonypy as sp
import yaml
from threadpoolctl import threadpool_limits


def write_fragment(frame, path):
    fragment = frame.copy()
    fragment.columns = fragment.columns.map(str)
    if "__index__" in fragment.columns or not fragment.columns.is_unique:
        raise ValueError("Fragment columns must be unique and must not contain reserved column '__index__'.")
    if (
        fragment.index.nlevels != 1
        or not fragment.index.is_unique
        or fragment.index.hasnans
        or not all(isinstance(value, str) for value in fragment.index)
    ):
        raise ValueError("Fragment axis identifiers must be unique, non-null strings.")
    fragment = fragment.reset_index(names="__index__")
    fragment.to_parquet(str(Path(path).resolve()), engine="pyarrow", index=False, compression="snappy")


threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
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

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")
write_fragment(pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names), f"X_{prefix}.parquet")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "symphonypy": importlib.metadata.version("symphonypy"),
        "pandas": pd.__version__,
        "pyarrow": pa.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
