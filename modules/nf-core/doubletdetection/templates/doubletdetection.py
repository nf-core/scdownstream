#!/usr/bin/env python3

import os
import platform
from pathlib import Path

os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ["NUMBA_CACHE_DIR"] = "./tmp"

import anndata as ad
import doubletdetection
import pyarrow
import yaml

adata = ad.read_h5ad("${h5ad}")

clf = doubletdetection.BoostClassifier(n_jobs=int("${task.cpus}"))
doublets = clf.fit(adata.X).predict()
scores = clf.doublet_score()

adata.obs["doublet"] = [label == 1 for label in doublets]
adata.obs["doublet_score"] = scores

adata.write_h5ad("${prefix}.h5ad", compression="gzip")

df = adata.obs[["doublet"]].copy()
df.columns = ["${prefix}"]
if df.index.hasnans or not df.index.is_unique or not all(isinstance(cell, str) for cell in df.index):
    raise ValueError("Doublet predictions require unique, non-null string cell IDs")
df.insert(0, "__index__", df.index)
df.to_parquet(Path("${prefix}.parquet").resolve(), index=False, engine="pyarrow", compression="snappy")
scores = adata.obs[["doublet_score"]].rename(columns={"doublet_score": "doubletdetection_score"})
scores.insert(0, "__index__", scores.index)
scores.to_parquet(Path("${prefix}.scores.parquet").resolve(), index=False, engine="pyarrow", compression="snappy")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "doubletdetection": doubletdetection.__version__,
        "pyarrow": pyarrow.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
