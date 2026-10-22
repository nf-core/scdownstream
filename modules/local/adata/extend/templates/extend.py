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


def restore_dtype(values, dtypes):
    categoricals = [dtype for dtype in dtypes if isinstance(dtype, pd.CategoricalDtype)]
    if categoricals:
        categories = categoricals[0].categories
        for dtype in categoricals[1:]:
            categories = categories.union(dtype.categories, sort=False)
        categories = categories.union(pd.Index(values.dropna().unique().tolist()), sort=False)
        return pd.Categorical(values, categories=categories, ordered=categoricals[0].ordered)

    values = values.infer_objects(copy=False)
    # AnnData cannot serialise object arrays containing booleans and missing values.
    if values.dtype == object and all(isinstance(value, (bool, np.bool_)) for value in values.dropna()):
        return pd.Categorical(values, categories=pd.Index([False, True], dtype=bool))
    return values


def merge_dataframe(base, paths):
    """Apply indexed updates while checking agreement between fragment sources."""
    fragments = []
    for path in paths:
        fragment = load_pickle_or_csv(path)
        if not fragment.index.is_unique:
            raise ValueError(f"Fragment '{path.name}' contains duplicate indices")
        if not fragment.columns.is_unique:
            raise ValueError(f"Fragment '{path.name}' contains duplicate column names")
        fragments.append((path, fragment.loc[fragment.index.isin(base.index)]))

    result = base.copy()
    columns = dict.fromkeys(column for _, fragment in fragments for column in fragment.columns)
    for column in columns:
        parts = {path.name: fragment[column] for path, fragment in fragments if column in fragment}
        incoming = pd.concat({name: values.dropna().astype(object) for name, values in parts.items()})
        grouped = incoming.groupby(level=1, sort=False, observed=True)
        conflicts = grouped.nunique().gt(1)
        if conflicts.any():
            rows = incoming.index.get_level_values(1).isin(conflicts.index[conflicts])
            sources = sorted(incoming.index.get_level_values(0)[rows].unique())
            raise ValueError(f"Conflicting non-null values in column {column!r} between fragments {sources}")

        values = (
            base[column].astype(object).copy() if column in base else pd.Series(np.nan, index=base.index, dtype=object)
        )
        updates = grouped.first()
        values.loc[updates.index] = updates
        values = values.where(values.notna(), np.nan)
        dtypes = ([base[column].dtype] if column in base else []) + [part.dtype for part in parts.values()]
        result[column] = restore_dtype(values, dtypes)
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
