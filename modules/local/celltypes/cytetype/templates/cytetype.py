#!/usr/bin/env python3
"""Per-sample CyteType annotation after lightweight Scanpy clustering + DE."""

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from importlib.metadata import PackageNotFoundError, version as pkg_version
from threadpoolctl import threadpool_limits
from cytetype import CyteType

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")
sc.settings.seed = 42
np.random.seed(42)

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
study_context = """${study_context}"""
if not study_context.strip():
    raise ValueError("cytetype_study_context must be a non-empty string when CyteType is enabled.")

symbol_col = "${symbol_col}"
leiden_resolution = float("${leiden_resolution}")
n_top_genes = int("${n_top_genes}")
auth_raw = """${auth_token}"""
auth_token_arg = auth_raw.strip() or None

orig_obs_cols = set(adata.obs.columns)
working = adata.copy()

if symbol_col != "index" and symbol_col:
    if symbol_col not in working.var.columns:
        raise ValueError(f"Symbol column {symbol_col} not found in adata.var.columns")
    working.var_names = working.var[symbol_col]

working.var_names = working.var_names.astype(str)

# QC-stage AnnData is typically raw counts in X — normalize + log for clustering / markers
sc.pp.normalize_total(working, target_sum=1e4)
sc.pp.log1p(working)

n_var = working.n_vars
n_top_hvg = max(500, min(2000, n_var))
sc.pp.highly_variable_genes(
    working,
    n_top_genes=n_top_hvg,
    flavor="seurat",
    subset=True,
)
sc.pp.scale(working, max_value=10)

n_obs = working.n_obs
n_pcs = min(50, max(1, n_obs - 1), max(1, working.n_vars - 1))
sc.tl.pca(working, n_comps=n_pcs)
n_neighbors = min(15, max(2, n_obs - 1))
sc.pp.neighbors(working, n_neighbors=n_neighbors, random_state=42)
sc.tl.umap(working, random_state=42)

group_key = "cytetype_leiden"
sc.tl.leiden(
    working,
    resolution=leiden_resolution,
    key_added=group_key,
    random_state=42,
)
sc.tl.rank_genes_groups(working, group_key, method="wilcoxon")

annotator = CyteType(
    working,
    group_key=group_key,
    rank_key="rank_genes_groups",
    n_top_genes=n_top_genes,
    auth_token=auth_token_arg,
)
working = annotator.run(
    study_context=study_context,
    show_progress=False,
    timeout_seconds=7200,
    require_artifacts=True,
)

added_cols = [c for c in working.obs.columns if c not in orig_obs_cols]
if not added_cols:
    raise RuntimeError(
        "CyteType did not add any new obs columns; check logs and study_context."
    )

df_out = working.obs[added_cols].reindex(adata.obs.index)
df_out.to_pickle(f"{prefix}.pkl")

adata.obs = pd.concat([adata.obs, df_out], axis=1)
adata.write_h5ad(f"{prefix}.h5ad")

try:
    cytetype_ver = pkg_version("cytetype")
except PackageNotFoundError:
    cytetype_ver = "unknown"

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "cytetype": cytetype_ver,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f, default_flow_style=False, sort_keys=True)
