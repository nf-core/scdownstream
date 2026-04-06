#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform
import anndata as ad
import scanpy as sc
from scipy.sparse import csc_matrix
import numpy as np
import scipy as sp
import yaml

adata = ad.read_h5ad("${h5ad}")

integration_methods = ["harmony", "scvi", "scanvi", "scimilarity", "seurat", "bbknn", "combat"]

for integration in integration_methods:
    embedding_key = f"X_{integration}"
    if embedding_key in adata.obsm.keys():
        adata.obsm[integration] = adata.obsm.pop(embedding_key)

for layer in adata.layers.keys():
    adata.layers[layer] = csc_matrix(adata.layers[layer]).astype(np.float32)
adata.X = csc_matrix(adata.X).astype(np.float32)
sc.pp.log1p(adata)

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scipy": sp.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
