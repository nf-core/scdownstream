#!/usr/bin/env python3

import platform
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow
import scanpy as sc
import yaml
from scipy.stats import entropy

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
cluster_col = "${cluster_col}"
annotation_col = "${annotation_col}"
integration = "${integration}"
resolution = "${resolution}"
subset = "${meta.subset ?: ''}"

if ":per_cell" not in annotation_col:
    raise ValueError(f"Unexpected per-cell annotation column: {annotation_col!r}")

output_col = annotation_col.replace(":per_cell", f":aggregated:{integration}:{resolution}")
if subset and subset != "global":
    output_col += f":subset:{subset}"

if annotation_col not in adata.obs.columns:
    raise ValueError(f"Annotation column {annotation_col!r} not found in adata.obs")

if cluster_col not in adata.obs.columns:
    raise ValueError(f"Cluster column {cluster_col!r} not found in adata.obs")

summary_rows = []
for cluster, group in adata.obs.groupby(cluster_col, observed=True):
    counts = group[annotation_col].astype(str).value_counts(normalize=True)
    dominant_label = counts.idxmax() if not counts.empty else ""
    dominant_fraction = float(counts.max()) if not counts.empty else 0.0
    label_entropy = float(entropy(counts, base=2)) if not counts.empty else 0.0
    conf_col = f"{annotation_col}:confidence"
    mean_conf = float(group[conf_col].mean()) if conf_col in group.columns else np.nan
    summary_rows.append(
        {
            cluster_col: cluster,
            "annotation_column": annotation_col,
            "output_column": output_col,
            "majority_label": dominant_label,
            "majority_fraction": dominant_fraction,
            "label_entropy": label_entropy,
            "mean_confidence": mean_conf,
            "n_cells": int(group.shape[0]),
        }
    )

cluster_majority = adata.obs.groupby(cluster_col, observed=True)[annotation_col].apply(
    lambda series: series.astype(str).value_counts().idxmax() if not series.empty else ""
)
obs_cols = {
    output_col: adata.obs[cluster_col].map(cluster_majority),
}

if summary_rows:
    pd.DataFrame(summary_rows).to_csv(f"{prefix}_aggregate_per_cell_annotations.csv", index=False)

obs_fragment = pd.DataFrame(obs_cols, index=adata.obs.index).copy()
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
adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pyarrow": pyarrow.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
