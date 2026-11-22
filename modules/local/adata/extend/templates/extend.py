#!/usr/bin/env python3

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


def load_parquet(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path, engine="pyarrow").set_index("__index__").rename_axis(None)


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

    if all(pd.api.types.is_bool_dtype(dtype) for dtype in dtypes):
        dtype = (
            "boolean" if values.isna().any() or any(isinstance(dtype, pd.BooleanDtype) for dtype in dtypes) else bool
        )
        return values.astype(dtype)

    if all(pd.api.types.is_numeric_dtype(dtype) for dtype in dtypes):
        dtype = np.result_type(*(getattr(dtype, "numpy_dtype", dtype) for dtype in dtypes))
        if all(pd.api.types.is_integer_dtype(source) for source in dtypes):
            if dtype.kind == "f":
                # Mixed signed and unsigned 64-bit integers have no exact NumPy common dtype.
                dtype = np.dtype("uint64" if values.dropna().ge(0).all() else "int64")
            nullable = values.isna().any() or any(
                isinstance(source, pd.api.extensions.ExtensionDtype) for source in dtypes
            )
            if nullable:
                dtype = f"{'UInt' if dtype.kind == 'u' else 'Int'}{dtype.itemsize * 8}"
            return values.astype(dtype)
        result = values.astype(dtype)
        if dtype.kind == "f" and any(pd.api.types.is_integer_dtype(source) for source in dtypes):
            for original, converted in zip(values, result, strict=True):
                if isinstance(original, (int, np.integer)) and not isinstance(original, (bool, np.bool_)):
                    if not np.isfinite(converted) or int(converted) != original:
                        raise ValueError("Merging numeric fragment dtypes would lose integer precision")
        return result

    if values.notna().any() and all(isinstance(value, (bool, np.bool_)) for value in values.dropna()):
        return values.astype("boolean")
    values = values.infer_objects(copy=False)
    return values


def merge_dataframe(base, paths):
    """Apply indexed updates while checking agreement between fragment sources."""
    fragments = []
    for path in paths:
        fragment = load_parquet(path)
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
    df = load_parquet(path).reindex(adata.obs_names)
    adata.obsm[path.stem] = np.float32(df.to_numpy())

for path in obsp_paths:
    adata.obsp[path.stem] = np.load(path, allow_pickle=True).item()

for path in uns_paths:
    if path.suffix == ".parquet":
        adata.uns[path.stem] = load_parquet(path)
    elif path.suffix == ".pkl":
        with path.open("rb") as handle:
            adata.uns[path.stem] = pickle.load(handle)
    else:
        raise ValueError(f"Unsupported uns file extension: {path}")

for path in layers_paths:
    layer = load_layer(path)
    if not issparse(layer):
        layer = np.asarray(layer, dtype=np.float32)
    adata.layers[path.stem] = layer

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")
adata.obs.to_csv(f"{prefix}_metadata.csv")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "pandas": pd.__version__,
        "pyarrow": importlib.metadata.version("pyarrow"),
        "numpy": np.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
