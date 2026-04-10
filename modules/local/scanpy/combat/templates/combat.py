#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import argparse
import shlex
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

args = argparse.Namespace(decimals=None)
if "${args}" != "":
    parser = argparse.ArgumentParser()
    parser.add_argument('--decimals', type=int, default=None)
    args = parser.parse_args(shlex.split("${args}"))

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

sc.pp.combat(adata, key="${batch_col}")
adata.X = csr_matrix(adata.X)

sc.pp.pca(adata)

if args.decimals is not None:
    adata.X.data = np.round(adata.X.data.astype(np.float64), args.decimals)
    adata.X.eliminate_zeros()
    adata.obsm["X_pca"] = np.round(adata.obsm["X_pca"].astype(np.float64), args.decimals)
    adata.varm["PCs"] = np.round(adata.varm["PCs"].astype(np.float64), args.decimals)
    adata.uns["pca"]["variance"] = np.round(adata.uns["pca"]["variance"].astype(np.float64), args.decimals)

adata.obsm["X_emb"] = adata.obsm["X_pca"]

adata.write_h5ad(f"{prefix}.h5ad")

np.save(f"{prefix}.npy", adata.X)

df = pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names)
df.to_pickle("X_${prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
