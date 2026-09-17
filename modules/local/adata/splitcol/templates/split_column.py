#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import importlib.metadata
import platform

import anndata as ad
import pandas as pd
import yaml


def _dataframe_for_h5ad(df: pd.DataFrame) -> pd.DataFrame:
    """Rewrite string indexes and columns so downstream readers can read the H5AD."""
    df = df.copy()
    if isinstance(df.index.dtype, pd.StringDtype):
        df.index = pd.CategoricalIndex(df.index.astype(str).to_list(), name=df.index.name)
    for col in df.columns:
        if isinstance(df[col].dtype, pd.StringDtype):
            df[col] = pd.Series(df[col].astype(str).to_list(), index=df.index, name=col, dtype=object)
    return df


def _prepare_adata_for_h5ad(adata_obj):
    """Avoid nullable-string-array encodings that downstream tools cannot read."""
    adata_obj.obs = _dataframe_for_h5ad(adata_obj.obs)
    adata_obj.var = _dataframe_for_h5ad(adata_obj.var)
    for uns_key, uns_value in list(adata_obj.uns.items()):
        if isinstance(uns_value, pd.DataFrame):
            adata_obj.uns[uns_key] = _dataframe_for_h5ad(uns_value)
        elif isinstance(uns_value, dict):
            for key, value in list(uns_value.items()):
                if isinstance(value, pd.DataFrame):
                    uns_value[key] = _dataframe_for_h5ad(value)


adata = ad.read_h5ad("${h5ad}")
column = "${column}"

assert column in adata.obs.columns, f"Column {column} not found in adata."

for value in adata.obs[column].unique():
    adata_subset = adata[adata.obs[column] == value].copy()
    value = value.replace(" ", "_")
    _prepare_adata_for_h5ad(adata_subset)
    adata_subset.write_h5ad(f"{value}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
