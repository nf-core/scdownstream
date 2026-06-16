#!/usr/bin/env python3

import os

os.environ["KMP_AFFINITY"] = "disabled"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"

import json
import math
import platform
import warnings
from dataclasses import replace

import numpy as np
import scanpy as sc
import scib_metrics
import yaml
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

prefix = "${prefix}"
h5ad_path = "${h5ad}"
max_cells_raw = "${task.ext.max_cells}"
subsample_strategy = "${task.ext.subsample_strategy}"
subsample_seed = int("${task.ext.subsample_seed}")
metric_profile = "${task.ext.metric_profile}"
neighbor_backend = "${task.ext.neighbor_backend}"
n_jobs = int("${task.cpus}")

max_cells = None
if max_cells_raw and max_cells_raw not in ("null", "None", ""):
    max_cells = int(max_cells_raw)


def _stratified_subsample(adata, n_max, strategy, seed, label_key="label", batch_key="batch"):
    n_before = adata.n_obs
    info = {
        "subsampled": False,
        "n_cells_before": n_before,
        "n_cells_after": n_before,
        "strategy": "none",
        "seed": seed,
        "max_cells": n_max,
    }
    if n_max is None or n_max <= 0 or n_before <= n_max:
        return adata, info

    if strategy in ("none", ""):
        strategy = "stratified_label_batch"
        warnings.warn(
            "scib_max_cells is set but scib_subsample_strategy is 'none'; "
            "using stratified_label_batch instead of uniform sampling."
        )

    rng = np.random.default_rng(seed)
    if strategy == "stratified_label":
        group_keys = adata.obs[label_key].astype(str)
    elif strategy == "stratified_label_batch":
        group_keys = (
            adata.obs[label_key].astype(str)
            + "\0"
            + adata.obs[batch_key].astype(str)
        )
    else:
        raise SystemExit(
            f"Unknown scib_subsample_strategy '{strategy}'; "
            "expected stratified_label or stratified_label_batch."
        )

    frac = n_max / n_before
    selected = []
    for group in group_keys.unique():
        idx = np.flatnonzero((group_keys == group).to_numpy())
        n_take = max(1, min(len(idx), int(round(len(idx) * frac))))
        if n_take >= len(idx):
            selected.extend(idx.tolist())
        else:
            selected.extend(rng.choice(idx, size=n_take, replace=False).tolist())

    if len(selected) > n_max:
        selected = rng.choice(np.array(selected), size=n_max, replace=False).tolist()

    adata_sub = adata[sorted(selected)].copy()
    info.update(
        {
            "subsampled": True,
            "n_cells_after": adata_sub.n_obs,
            "strategy": strategy,
        }
    )
    return adata_sub, info


def _metric_config(profile):
    if profile == "full":
        return None, None
    if profile != "fast":
        raise SystemExit(
            f"Unknown scib_metric_profile '{profile}'; expected 'fast' or 'full'."
        )
    return (
        replace(
            BioConservation(),
            isolated_labels=False,
            nmi_ari_cluster_labels_kmeans=False,
        ),
        replace(BatchCorrection(), pcr_comparison=False),
    )


def _neighbor_computer(backend):
    if backend in ("default", ""):
        return None

    if backend != "faiss":
        raise SystemExit(
            f"Unknown scib_neighbor_backend '{backend}'; expected 'default' or 'faiss'."
        )

    import faiss
    from scib_metrics.nearest_neighbors import NeighborsResults

    def faiss_brute_force_nn(X: np.ndarray, k: int):
        X = np.ascontiguousarray(X, dtype=np.float32)
        index = faiss.IndexFlatL2(X.shape[1])
        index.add(X)
        distances, indices = index.search(X, k)
        return NeighborsResults(indices=indices, distances=np.sqrt(distances))

    return faiss_brute_force_nn


adata = sc.read_h5ad(h5ad_path)

missing = [c for c in ("batch", "label") if c not in adata.obs]
if missing:
    raise SystemExit(
        f"scib-metrics requires obs columns {missing}; available: {list(adata.obs.columns)}"
    )
if "X_emb" not in adata.obsm:
    raise SystemExit(
        f"scib-metrics requires obsm['X_emb']; available: {list(adata.obsm.keys())}"
    )

if "highly_variable" not in adata.var.columns:
    adata.var["highly_variable"] = True
elif not bool(adata.var["highly_variable"].any()):
    warnings.warn(
        "No highly_variable genes flagged; using all genes as HVG for scib-metrics PCA."
    )
    adata.var["highly_variable"] = True

adata, subsample_info = _stratified_subsample(
    adata,
    max_cells,
    subsample_strategy,
    subsample_seed,
)

labels = adata.obs["label"].astype(str)
if (labels == "Unknown").all():
    warnings.warn(
        "All cells have label 'Unknown'; bio-conservation metrics are not meaningful."
    )

ad_bm = adata.copy()
bio_metrics, batch_metrics = _metric_config(metric_profile)
skip_pcr = metric_profile == "fast"

if not skip_pcr:
    if "counts" not in ad_bm.layers:
        ad_bm.layers["counts"] = ad_bm.X.copy()
    ad_bm.X = ad_bm.layers["counts"].copy()
    sc.pp.normalize_total(ad_bm, target_sum=1e4)
    sc.pp.log1p(ad_bm)
elif ad_bm.X.max() > 30:
    warnings.warn(
        "scib fast profile skips PCR comparison; assuming adata.X is already normalized."
    )

bm_kw = {"n_jobs": n_jobs}
if bio_metrics is not None:
    bm_kw["bio_conservation_metrics"] = bio_metrics
if batch_metrics is not None:
    bm_kw["batch_correction_metrics"] = batch_metrics

if labels.nunique() <= 1:
    warnings.warn(
        "obs['label'] has only one unique value; disabling BRAS (silhouette batch)."
    )
    batch_cfg = bm_kw.get("batch_correction_metrics", BatchCorrection())
    bm_kw["batch_correction_metrics"] = replace(batch_cfg, bras=False)

neighbor_computer = _neighbor_computer(neighbor_backend)

bm = Benchmarker(
    ad_bm,
    batch_key="batch",
    label_key="label",
    embedding_obsm_keys=["X_emb"],
    pre_integrated_embedding_obsm_key=None,
    progress_bar=False,
    **bm_kw,
)
bm.prepare(neighbor_computer=neighbor_computer)
bm.benchmark()
results = bm.get_results(min_max_scale=False, clean_names=True)
results.to_csv(f"{prefix}_metrics.tsv", sep="\t")

benchmark_info = {
    "integration_method": prefix,
    "metric_profile": metric_profile,
    "neighbor_backend": neighbor_backend or "default",
    "n_jobs": n_jobs,
    **subsample_info,
}
with open(f"{prefix}_benchmark_info.json", "w") as f_info:
    json.dump(benchmark_info, f_info, indent=2)


def _mqc_table_cell(v):
    if v is None:
        return None
    try:
        fv = float(v)
        if math.isnan(fv) or math.isinf(fv):
            return None
        return fv
    except (TypeError, ValueError):
        pass
    if hasattr(v, "item"):
        return _mqc_table_cell(v.item())
    return v


def _benchmark_description(info):
    parts = [
        "scib-metrics benchmark (one row per integration method).",
        f"Profile: {info['metric_profile']}.",
        f"Neighbor backend: {info['neighbor_backend']}.",
    ]
    if info["subsampled"]:
        parts.append(
            "Subsampled "
            f"{info['n_cells_before']} -> {info['n_cells_after']} cells "
            f"({info['strategy']}, seed={info['seed']})."
        )
    else:
        parts.append(f"Cells: {info['n_cells_before']} (no subsampling).")
    parts.append(
        "Subsampled scores are for within-pipeline monitoring; "
        "not directly comparable to full-data scIB paper benchmarks."
    )
    return " ".join(parts)


# MultiQC custom table: data is { sample (row): { column: value, ... } }.
# Fixed section id merges all integration runs (one row each); one embedding row per run.
metrics_row = results.iloc[0]
row = {str(c): _mqc_table_cell(metrics_row[c]) for c in results.columns}
mqc_data = {prefix: row}
mqc_headers = {str(c): {"format": "{:.3f}"} for c in results.columns}

with open(f"{prefix}_mqc.json", "w") as f_json:
    json.dump(
        {
            "id": "scib_metrics_benchmark",
            "plot_type": "table",
            "section_name": "scib-metrics",
            "description": _benchmark_description(benchmark_info),
            "pconfig": {"col1_header": "Integration method"},
            "headers": mqc_headers,
            "data": mqc_data,
        },
        f_json,
    )


versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scib-metrics": scib_metrics.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
