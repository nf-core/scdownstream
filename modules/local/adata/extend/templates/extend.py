#!/usr/bin/env python3

import platform
import os
import pickle

import anndata as ad
import pandas as pd
import numpy as np
import yaml

adata = ad.read_h5ad("${base}")
prefix = "${prefix}"
obs_paths = sorted("${obs}".split())
var_paths = sorted("${var}".split())
obsm_paths = sorted("${obsm}".split())
obsp_paths = sorted("${obsp}".split())
uns_paths = sorted("${uns}".split())
layers_paths = sorted("${layers}".split())

def simple_name(path):
    basename = os.path.basename(path)
    return basename[:basename.rfind(".")]

def load_pickle_or_csv(path):
    if path.endswith(".pkl"):
        return pd.read_pickle(path)
    elif path.endswith(".csv"):
        return pd.read_csv(path, index_col = 0)
    else:
        raise ValueError(f"Unsupported file extension: {path}")

for path in obs_paths:
    df = load_pickle_or_csv(path).reindex(adata.obs_names)
    adata.obs = pd.concat([adata.obs, df], axis=1)

for path in var_paths:
    df = load_pickle_or_csv(path).reindex(adata.var_names)
    adata.var = pd.concat([adata.var, df], axis=1)

# Ensure 'highly_variable' is boolean if present
if "highly_variable" in adata.var:
    adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)

for path in obsm_paths:
    df = pd.read_pickle(path).reindex(adata.obs_names)
    adata.obsm[simple_name(path)] = np.float32(df.to_numpy())

for path in obsp_paths:
    adata.obsp[simple_name(path)] = np.load(path, allow_pickle=True).item()

for path in uns_paths:
    adata.uns[simple_name(path)] = pickle.load(open(path, "rb"))

for path in layers_paths:
    adata.layers[simple_name(path)] = np.float32(np.load(path))

adata.write_h5ad(f"{prefix}.h5ad")
adata.obs.to_csv(f"{prefix}_metadata.csv")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
