#!/usr/bin/env python3

import os

# These are needed to prevent errors during import of scanpy
# when using singularity/apptainer
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform
from pathlib import Path

import pyarrow
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
batch_col = "${batch_col ?: ''}"

kwargs = {}
if batch_col and adata.obs[batch_col].nunique() > 1:
    kwargs["batch_key"] = batch_col

sc.pp.scrublet(adata, **kwargs)

adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
df = adata.obs[["predicted_doublet"]].copy()
df.columns = ["${prefix}"]
if df.index.hasnans or not df.index.is_unique or not all(isinstance(cell, str) for cell in df.index):
    raise ValueError("Doublet predictions require unique, non-null string cell IDs")
df.insert(0, "__index__", df.index)
df.to_parquet(Path("${prefix}.parquet").resolve(), index=False, engine="pyarrow", compression="snappy")

adata = adata[~adata.obs["predicted_doublet"]].copy()

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pyarrow": pyarrow.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
