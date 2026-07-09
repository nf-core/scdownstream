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

adata_proc = adata.copy()
sc.experimental.pp.normalize_pearson_residuals(adata_proc, layer="counts")

layer_matrix = adata_proc.X.copy()
adata.layers["pearson_residuals"] = layer_matrix

if issparse(layer_matrix):
    save_npz("pearson_residuals.npz", layer_matrix.astype("float32"))
else:
    np.save("pearson_residuals.npy", np.asarray(layer_matrix, dtype=np.float32))

adata.write_h5ad(f"{prefix}.h5ad")

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
