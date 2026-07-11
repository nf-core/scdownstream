#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"

import numpy as np
import scanpy as sc
import yaml
from scipy.sparse import issparse, save_npz

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

adata_norm = adata.copy()
adata_norm.X = adata.layers["counts"]
sc.pp.normalize_total(adata_norm, target_sum=None)
sc.pp.log1p(adata_norm)

layer_matrix = adata_norm.X.copy()
adata.layers["log1p"] = layer_matrix

if issparse(layer_matrix):
    save_npz("log1p.npz", layer_matrix.astype("float32"))
else:
    np.save("log1p.npy", np.asarray(layer_matrix, dtype=np.float32))

adata.write_h5ad(f"{prefix}.h5ad")

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
