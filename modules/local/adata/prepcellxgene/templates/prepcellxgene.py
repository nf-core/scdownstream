#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform

import anndata as ad
import numpy as np
import scanpy as sc
import scipy as sp
import yaml
from scipy.sparse import csc_matrix

adata = ad.read_h5ad("${h5ad}")

integration_methods = ["scanvi", "scvi", "symphony", "scimilarity", "seurat", "bbknn", "combat", "pca", "expimap"]

dim_reds = ["umap", "tsne"]

for key in list(adata.obsm.keys()):
    if not any(key.endswith(dim_red) for dim_red in dim_reds):
        del adata.obsm[key]

# Delete everything in uns and layers (CELLxGENE can't display these anyway)
adata.uns = {}
adata.layers = {}

# Convert all float64 columns to float32
for df in [adata.obs, adata.var]:
    for col in df.columns:
        if df[col].dtype == np.float64:
            df[col] = df[col].astype(np.float32)

adata.X = csc_matrix(adata.X).astype(np.float32)
sc.pp.log1p(adata)

# Merged uns may still contain pandas StringDtype indexes (e.g. scanpy pts).
# Prefer source-side CategoricalIndex sanitisation; keep True as a write safety net.
ad.settings.allow_write_nullable_strings = True
adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scipy": sp.__version__,
        "numpy": np.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
