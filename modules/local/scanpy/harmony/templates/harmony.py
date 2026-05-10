#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import yaml
import argparse
import shlex

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import harmonypy
import scanpy as sc
import pandas as pd
import numpy as np

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
args = "${args}"
parser = argparse.ArgumentParser()
parser.add_argument("--decimals", type=int, default=None)
params = parser.parse_args(shlex.split(args))

prefix = "${prefix}"

adata_processing = adata.copy()

if "${counts_layer}" != "X":
    adata_processing.X = adata.layers["${counts_layer}"]

sc.pp.log1p(adata_processing)
sc.pp.pca(adata_processing)

harmony_out = harmonypy.run_harmony(
    adata_processing.obsm["X_pca"].astype("float64"),
    adata_processing.obs,
    "${batch_col}",
)

emb = harmony_out.Z_corr

# harmonypy 0.2.0 changed Z_corr orientation; accept either layout.
# See https://github.com/potulabe/symphonypy/issues/8
if emb.shape == adata_processing.obsm["X_pca"].shape:
    adata_processing.obsm["X_emb"] = emb
elif emb.T.shape == adata_processing.obsm["X_pca"].shape:
    adata_processing.obsm["X_emb"] = emb.T
else:
    raise ValueError(
        f"Unexpected Harmony embedding shape {emb.shape}; "
        f"expected {adata_processing.obsm['X_pca'].shape} or its transpose."
    )

if params.decimals is not None:
    adata_processing.obsm["X_emb"] = adata_processing.obsm["X_emb"].round(params.decimals)
adata.obsm["X_emb"] = adata_processing.obsm["X_emb"]

adata.write_h5ad(f"{prefix}.h5ad")

df = pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "harmonypy": harmonypy.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
