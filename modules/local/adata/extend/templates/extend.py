#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import pickle
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.io import mmread
from scipy.sparse import issparse, load_npz

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
        return pd.read_csv(path, index_col=0)
    else:
        raise ValueError(f"Unsupported file extension: {path}")


def load_layer(path):
    if path.suffix == ".npz":
        return load_npz(path)
    if path.suffix == ".mtx":
        return mmread(path).tocsr().astype(np.float32)
    if path.suffix == ".npy":
        return np.float32(np.load(path))
    raise ValueError(f"Unsupported layer file extension: {path}")


def _dataframe_without_nullable_strings(df: pd.DataFrame) -> pd.DataFrame:
    """Cast pandas StringDtype columns and indexes to plain object strings."""
    df = df.copy()
    if isinstance(df.index.dtype, pd.StringDtype):
        df.index = pd.Index(df.index.to_numpy(dtype=object), name=df.index.name)
    for col in df.columns:
        if isinstance(df[col].dtype, pd.StringDtype):
            df[col] = df[col].astype(object)
    return df


def _prepare_adata_for_h5ad(adata_obj):
    adata_obj.obs = _dataframe_without_nullable_strings(adata_obj.obs)
    adata_obj.var = _dataframe_without_nullable_strings(adata_obj.var)
    for uns_key, uns_value in list(adata_obj.uns.items()):
        if isinstance(uns_value, pd.DataFrame):
            adata_obj.uns[uns_key] = _dataframe_without_nullable_strings(uns_value)
        elif isinstance(uns_value, dict):
            for key, value in list(uns_value.items()):
                if isinstance(value, pd.DataFrame):
                    uns_value[key] = _dataframe_without_nullable_strings(value)


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
    layer = load_layer(path)
    if not issparse(layer):
        layer = np.asarray(layer, dtype=np.float32)
    adata.layers[path.stem] = layer

_prepare_adata_for_h5ad(adata)
adata.write_h5ad(f"{prefix}.h5ad")
adata.obs.to_csv(f"{prefix}_metadata.csv")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "pandas": pd.__version__,
        "numpy": np.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
