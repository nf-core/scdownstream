#!/usr/bin/env python3
"""Per-sample CyteType annotation on pre-computed clustering and DE."""

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from importlib.metadata import version
from cytetype import CyteType

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
study_context = "${study_context}"
if not study_context.strip():
    raise ValueError("cytetype_study_context must be a non-empty string when CyteType is enabled.")

symbol_col = "${symbol_col}"
group_key = "${group_key}"
rank_key = "${rank_key}"
_auth = os.environ.get("CYTETYPE_API_KEY")
auth_token_arg = _auth.strip() if _auth and _auth.strip() else None

orig_obs_cols = set(adata.obs.columns)
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
adata_work = annotator.run(
    study_context=study_context,
    show_progress=False,
    timeout_seconds=7200,
    require_artifacts=True,
)

added_cols = [c for c in adata_work.obs.columns if c not in orig_obs_cols]
if not added_cols:
    raise RuntimeError(
        "CyteType did not add any new obs columns; check logs and study_context."
    )

df_out = adata_work.obs[added_cols].reindex(adata.obs.index)
df_out.to_pickle(f"{prefix}.pkl")

adata.obs = pd.concat([adata.obs, df_out], axis=1)
adata.write_h5ad(f"{prefix}.h5ad")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "cytetype": version("cytetype")
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f, default_flow_style=False, sort_keys=True)
