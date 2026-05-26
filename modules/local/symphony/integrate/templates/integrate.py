#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform
import yaml

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import symphonypy as sp
import pandas as pd

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")

prefix = "${prefix}"

adata_processing = adata.copy()

if "${counts_layer}" != "X":
    adata_processing.X = adata.layers["${counts_layer}"]

sc.pp.log1p(adata_processing)
sc.pp.pca(adata_processing)

sp.pp.harmony_integrate(
    adata_processing,
    key="${batch_col}",
    flavor="python",
    ref_basis_source="X_pca",
    ref_basis_adjusted="X_pca_symphony",
)

adata.obsm["X_pca_symphony"] = adata_processing.obsm["X_pca_symphony"]
adata.obsm["X_emb"] = adata_processing.obsm["X_pca_symphony"]
adata.uns["symphony"] = adata_processing.uns["harmony"]

adata.write_h5ad(f"{prefix}.h5ad")

df = pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "symphonypy": importlib.metadata.version("symphonypy"),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
