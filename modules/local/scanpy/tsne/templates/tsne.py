#!/usr/bin/env python3

import platform
from pathlib import Path

import pandas as pd
import pyarrow as pa
import scanpy as sc
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
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}", backed="r")
prefix = "${prefix}"
representation = "${representation}"

if not representation:
    raise ValueError("representation parameter is required but was not provided")

if representation not in adata.obsm:
    raise ValueError(f"Representation '{representation}' not found in adata.obsm. Available: {list(adata.obsm.keys())}")

sc.tl.tsne(adata, use_rep=representation, random_state=0)

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")
df = pd.DataFrame(adata.obsm["X_tsne"], index=adata.obs_names)
write_fragment(df, f"X_{prefix}.parquet")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__,
        "pyarrow": pa.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
