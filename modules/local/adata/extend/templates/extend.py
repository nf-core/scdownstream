#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import pickle
import importlib.metadata

import anndata as ad
import pandas as pd
import numpy as np
import yaml
from pathlib import Path

adata = ad.read_h5ad("${base}")
prefix = "${prefix}"
obs_paths = sorted(Path("obs/").glob("*"))
var_paths = sorted(Path("var/").glob("*"))
obsm_paths = sorted(Path("obsm/").glob("*"))
obsp_paths = sorted(Path("obsp/").glob("*"))
uns_paths = sorted(Path("uns/").glob("*"))
layers_paths = sorted(Path("layers/").glob("*"))

def load_pickle_or_csv(path):
    if path.suffix == ".pkl":
        return pd.read_pickle(path)
    elif path.suffix == ".csv":
        return pd.read_csv(path, index_col = 0)
    else:
        raise ValueError(f"Unsupported file extension: {path}")

for path in obs_paths:
    df = load_pickle_or_csv(path).reindex(adata.obs_names)
    adata.obs = pd.concat([adata.obs, df], axis=1)

for path in var_paths:
    df = load_pickle_or_csv(path).reindex(adata.var_names)
    adata.var = pd.concat([adata.var, df], axis=1)

for path in obsm_paths:
    df = pd.read_pickle(path).reindex(adata.obs_names)
    adata.obsm[path.stem] = np.float32(df.to_numpy())

for path in obsp_paths:
    adata.obsp[path.stem] = np.load(path, allow_pickle=True).item()

for path in uns_paths:
    adata.uns[path.stem] = pickle.load(open(path, "rb"))

for path in layers_paths:
    adata.layers[path.stem] = np.float32(np.load(path))

adata.write_h5ad(f"{prefix}.h5ad")
adata.obs.to_csv(f"{prefix}_metadata.csv")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "pandas": pd.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
