#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import pandas as pd
import scanpy as sc
import yaml
from scipy.sparse import csr_matrix
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
input_layer = "${input_layer}"
log_normalize = "${log_normalize}" == "true"
batch_col = "${batch_col}"

adata_proc = adata.copy()

if input_layer != "X" and input_layer not in adata_proc.layers:
    raise ValueError(
        f"input_layer {input_layer!r} is not present in adata.layers "
        f"(available: {list(adata_proc.layers.keys())})"
    )

if input_layer != "X":
    adata_proc.X = adata_proc.layers[input_layer]

if log_normalize:
    sc.pp.normalize_total(adata_proc)
    sc.pp.log1p(adata_proc)

combat_layer = csr_matrix(sc.pp.combat(adata_proc, key=batch_col, inplace=False))
adata_proc.layers["combat"] = combat_layer
sc.pp.pca(adata_proc, layer="combat")

adata.layers["combat"] = combat_layer
adata.obsm["X_emb"] = adata_proc.obsm["X_pca"]
adata.write_h5ad(f"{prefix}.h5ad")

pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names).to_pickle(f"X_{prefix}.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
