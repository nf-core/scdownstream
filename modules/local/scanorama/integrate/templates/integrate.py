#!/usr/bin/env python3

import os

os.environ["KMP_AFFINITY"] = "disabled"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import importlib.metadata
import platform

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
adata_proc = adata.copy()
prefix = "${prefix}"
batch_col = "${batch_col}"
input_layer = "${input_layer}"
log_normalize = "${log_normalize}" == "true"

if input_layer != "X" and input_layer not in adata_proc.layers:
    raise ValueError(
        f"input_layer {input_layer!r} is not present in adata.layers "
        f"(available: {list(adata_proc.layers.keys())})"
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
adata.write_h5ad(f"{prefix}.h5ad")

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
