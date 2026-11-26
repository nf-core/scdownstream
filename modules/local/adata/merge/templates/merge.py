#!/usr/bin/env python3

import os
import platform
from collections import defaultdict

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import yaml
from scipy.sparse import csr_matrix

input_paths = sorted("${h5ads}".split())
query_adatas = [sc.read_h5ad(f) for f in input_paths]
adatas = query_adatas

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

gene_presence_key = "gene_presence"
sample_ids = []
for adata in query_adatas:
    if adata.n_obs == 0:
        raise ValueError("Cannot merge an incoming AnnData with no cells.")
    sample_ids.append(str(adata.obs["sample"].iat[0]))

obs_col_intersection = set(adatas[0].obs.columns).intersection(*[adata.obs.columns for adata in adatas[1:]])
force_obs_cols = "${force_obs_cols}"
obs_col_intersection = sorted(obs_col_intersection.union(force_obs_cols.split(",") if force_obs_cols else []))


def get_columns(adata):
    return {dtype: adata.obs.select_dtypes(include=dtype).columns for dtype in ["object", "category", "number", "bool"]}


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
        if column_dtypes[col] in ({"number"}, {"bool"}):
            continue
        adata.obs[col] = adata.obs[col].astype(str).astype("category")

adata_outer = ad.concat(adatas, join="outer")
genes = adata_outer.var_names
presence = pd.DataFrame(
    np.column_stack([genes.isin(adata.var_names) for adata in query_adatas]),
    index=genes,
    columns=sample_ids,
)
presence_complete = not base_path
if base_path and gene_presence_key in adata_base.varm:
    base_presence = adata_base.varm[gene_presence_key].reindex(genes, fill_value=False)
    presence = pd.concat([base_presence, presence], axis=1)
    presence_complete = bool(adata_base.uns.get("gene_presence_complete", False))

if not presence.columns.is_unique:
    duplicates = presence.columns[presence.columns.duplicated()].unique().tolist()
    raise ValueError(f"Sample IDs must be unique across merge inputs. Repeated IDs: {duplicates}")

intersection = presence.all(axis=1).to_numpy()
if base_path:
    # Older atlases may preserve only the aggregate intersection.
    base_genes = adata_base.var_names
    if "intersection" in adata_base.var:
        base_genes = base_genes[adata_base.var["intersection"].eq(True).fillna(False)]
    intersection = intersection & genes.isin(base_genes)

if genes.name in presence.columns:
    # H5AD stores the feature index and sample columns under sibling keys.
    adata_outer.var_names = genes.rename(None)
adata_outer.varm[gene_presence_key] = presence
adata_outer.uns["gene_presence_complete"] = presence_complete
adata_outer.var["intersection"] = intersection
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
