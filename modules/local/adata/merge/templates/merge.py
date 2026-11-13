#!/usr/bin/env python3

import os
import platform
from collections import defaultdict

import anndata as ad
import numpy as np
import scanpy as sc
import scipy
import yaml
from scipy.sparse import csr_matrix

input_paths = sorted("${h5ads}".split())
adatas = [sc.read_h5ad(f) for f in input_paths]

base_path = "${base}"
if base_path:
    adata_base = sc.read_h5ad(base_path)
    adatas = [adata_base] + adatas
    input_paths = [base_path] + input_paths

provenance_columns = ["ambient_correction_method", "analysis_counts_source"]
correction_layers = {
    layer
    for adata in adatas
    for layer in adata.layers
    if isinstance(layer, str) and layer.startswith("ambient_corrected_")
}

for source, adata in zip(input_paths, adatas, strict=True):
    if "raw_counts" not in adata.layers:
        raise ValueError(f"AnnData '{source}' has no raw_counts layer.")
    adata.layers["raw_counts"] = csr_matrix(adata.layers["raw_counts"])
    for layer in sorted(correction_layers):
        if layer not in adata.layers:
            adata.layers[layer] = adata.layers["raw_counts"].copy()
        else:
            adata.layers[layer] = csr_matrix(adata.layers[layer])
    for column in provenance_columns:
        if column not in adata.obs:
            adata.obs[column] = "unknown"

genes = [adata.var_names for adata in adatas]
obs_col_intersection = set(adatas[0].obs.columns).intersection(*[adata.obs.columns for adata in adatas[1:]])
force_obs_cols = "${force_obs_cols}"
obs_col_intersection = sorted(obs_col_intersection.union(force_obs_cols.split(",") if force_obs_cols else []))


def get_columns(adata):
    return {dtype: adata.obs.select_dtypes(include=dtype).columns for dtype in ["object", "category", "number"]}


column_dtypes = defaultdict(set)

for adata in adatas:
    for dtype, columns in get_columns(adata).items():
        for column in columns:
            column_dtypes[column].add(dtype)

column_defaults = {}

for column, dtypes in column_dtypes.items():
    if len(dtypes) > 1:
        for adata in adatas:
            if column in adata.obs.columns:
                adata.obs[column] = adata.obs[column].astype("object")
        column_defaults[column] = "unknown"
    else:
        column_defaults[column] = np.nan if dtypes.copy().pop() == "number" else "unknown"

column_defaults["label"] = "Unknown"

for adata in adatas:
    for col in set(obs_col_intersection).difference(adata.obs.columns):
        adata.obs[col] = column_defaults[col]
    adata.obs = adata.obs[obs_col_intersection]

    for col in obs_col_intersection:
        if column_dtypes[col] == {"number"}:
            continue
        adata.obs[col] = adata.obs[col].astype(str).astype("category")

gene_intersection = set(genes[0]).intersection(*genes[1:])

adata_outer = ad.concat(adatas, join="outer")
adata_outer.var["intersection"] = adata_outer.var_names.isin(gene_intersection).astype(bool)
adata_outer.obs_names_make_unique()
adata_outer.X = csr_matrix(adata_outer.X)

# Sort obs columns alphabetically to make reproducible
adata_outer.obs = adata_outer.obs.reindex(sorted(adata_outer.obs.columns), axis=1)

adata_outer.write_h5ad("${prefix}_outer.h5ad", compression="gzip")

if base_path:
    adata_integrate = adata_outer[~adata_outer.obs.index.isin(adata_base.obs.index)].copy()

    known_labels = adata_base.obs["label"].unique()
    adata_integrate.obs["label"] = adata_integrate.obs["label"].map(lambda x: x if x in known_labels else "Unknown")
    adata_integrate.write_h5ad("${prefix}_integrate.h5ad", compression="gzip")
else:
    os.symlink("${prefix}_outer.h5ad", "${prefix}_integrate.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "scipy": scipy.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
