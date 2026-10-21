#!/usr/bin/env python3

import csv
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
        with path.open(newline="") as handle:
            columns = next(csv.reader(handle))[1:]
        if len(columns) != len(set(columns)):
            raise ValueError(f"Fragment '{path.name}' contains duplicate column names")
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


def merge_dataframe(base, paths):
    """Apply indexed updates while checking agreement between fragment sources."""
    updates = {}
    sources = {}
    categorical_dtypes = {}
    for path in paths:
        fragment = load_pickle_or_csv(path)
        if not fragment.index.is_unique:
            raise ValueError(f"Fragment '{path.name}' contains duplicate indices")
        if not fragment.columns.is_unique:
            raise ValueError(f"Fragment '{path.name}' contains duplicate column names")
        for column in fragment.columns:
            incoming = fragment[column].reindex(base.index)
            if column not in updates:
                updates[column] = pd.Series(np.nan, index=base.index, dtype=object)
                sources[column] = pd.Series(index=base.index, dtype=object)
                categorical_dtypes[column] = []
                if column in base and isinstance(base[column].dtype, pd.CategoricalDtype):
                    categorical_dtypes[column].append(base[column].dtype)
            if isinstance(incoming.dtype, pd.CategoricalDtype):
                categorical_dtypes[column].append(incoming.dtype)
            present = incoming.notna()
            overlap = present & updates[column].notna()
            conflicts = ~updates[column].loc[overlap].eq(incoming.loc[overlap].astype(object))
            if conflicts.any():
                previous = sorted(sources[column].loc[conflicts[conflicts].index].unique())
                raise ValueError(
                    f"Conflicting non-null values in column {column!r} between fragments {previous} and {path.name!r}"
                )
            updates[column].loc[present] = incoming.loc[present].astype(object)
            sources[column].loc[present] = path.name

    result = base.copy()
    for column, incoming in updates.items():
        values = (
            base[column].astype(object).copy() if column in base else pd.Series(np.nan, index=base.index, dtype=object)
        )
        present = incoming.notna()
        values.loc[present] = incoming.loc[present]
        values = values.where(values.notna(), np.nan)
        dtypes = categorical_dtypes[column]
        if dtypes:
            categories = dtypes[0].categories.copy()
            for dtype in dtypes[1:]:
                categories = categories.union(dtype.categories, sort=False)
            categories = categories.union(pd.Index(values.dropna().unique().tolist()), sort=False)
            result[column] = pd.Categorical(values, categories=categories, ordered=dtypes[0].ordered)
        else:
            values = values.infer_objects(copy=False)
            # AnnData cannot serialise object arrays containing booleans and missing values.
            if values.dtype == object and all(isinstance(value, (bool, np.bool_)) for value in values.dropna()):
                values = pd.Series(
                    pd.Categorical(values, categories=pd.Index([False, True], dtype=bool)), index=base.index
                )
            result[column] = values
    return result


adata.obs = merge_dataframe(adata.obs, obs_paths)
adata.var = merge_dataframe(adata.var, var_paths)

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
