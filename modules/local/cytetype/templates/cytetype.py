#!/usr/bin/env python3
"""Per-sample CyteType annotation on pre-computed clustering and DE."""

import os
import platform
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow
import scanpy as sc
import yaml
from cytetype import CyteType

adata = sc.read_h5ad("input.h5ad")
prefix = "${prefix}"
study_context = "${study_context}"
if not study_context.strip():
    raise ValueError("cytetype_study_context must be a non-empty string when CyteType is enabled.")

integration = "${integration}"
resolution = "${resolution}"
symbol_col = "${symbol_col}"
group_key = "${group_key}"
rank_key = "${rank_key}"
_auth = os.environ.get("CYTETYPE_API_KEY")
auth_token_arg = _auth.strip() if _auth and _auth.strip() else None

cytetype_base = f"annotation:cytetype:{integration}:{resolution}:{group_key}"
output_cols = [
    f"{cytetype_base}:cell_type",
    f"{cytetype_base}:ontology_term",
    f"{cytetype_base}:ontology_term_id",
    f"{cytetype_base}:cell_state",
]
library_cols = [
    f"cytetype_annotation_{group_key}",
    f"cytetype_cellOntologyTerm_{group_key}",
    f"cytetype_cellOntologyTermID_{group_key}",
    f"cytetype_cellState_{group_key}",
]

adata_work = adata.copy()

if symbol_col != "index" and symbol_col:
    if symbol_col not in adata_work.var.columns:
        raise ValueError(f"Symbol column {symbol_col} not found in adata.var.columns")
    adata_work.var_names = adata_work.var[symbol_col]

adata_work.var_names = adata_work.var_names.astype(str)

if group_key not in adata_work.obs.columns:
    raise ValueError(f"Group key {group_key!r} not found in adata.obs.columns")
if rank_key not in adata_work.uns:
    raise ValueError(f"Rank key {rank_key!r} not found in adata.uns")

annotator = CyteType(
    adata_work,
    group_key=group_key,
    rank_key=rank_key,
    auth_token=auth_token_arg,
)
# The API requires arrays even when the input has no visualisation coordinates.
if annotator.visualization_data.get("coordinates") is None:
    annotator.visualization_data = {"coordinates": [], "clusters": []}

adata_work = annotator.run(
    study_context=study_context,
    show_progress=False,
    timeout_seconds=7200,
    require_artifacts=True,
)

missing_cols = [c for c in library_cols if c not in adata_work.obs.columns]
if missing_cols:
    raise RuntimeError(f"CyteType did not add expected obs columns: {missing_cols}. Check logs and study_context.")

df_out = adata_work.obs[library_cols].rename(columns=dict(zip(library_cols, output_cols))).reindex(adata.obs.index)
obs_fragment = df_out.copy()
if "__index__" in obs_fragment.columns or not obs_fragment.columns.is_unique:
    raise ValueError("Observation fragment columns must be unique and must not contain reserved '__index__'")
if (
    obs_fragment.index.hasnans
    or not obs_fragment.index.is_unique
    or not all(isinstance(value, str) for value in obs_fragment.index)
):
    raise ValueError("Observation fragment index must contain unique, non-null strings")
obs_fragment.reset_index(names="__index__").to_parquet(
    Path(f"{prefix}.parquet").resolve(), engine="pyarrow", index=False, compression="snappy"
)

adata.obs = pd.concat([adata.obs, df_out], axis=1)
adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "pyarrow": pyarrow.__version__,
        "scanpy": sc.__version__,
        "cytetype": version("cytetype"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f, default_flow_style=False, sort_keys=True)
