#!/usr/bin/env python3

import os

os.environ["KMP_AFFINITY"] = "disabled"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"

import argparse
import json
import math
import platform
import shlex
import warnings
from dataclasses import replace

import scanpy as sc
import scib_metrics
import yaml
from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark._core import BatchCorrection

prefix = "${prefix}"
args_str = "${args}"
integration = "${integration_name}"
h5ad_path = "${h5ad}"

parser = argparse.ArgumentParser()
parser.add_argument("--n-jobs", type=int, default=1, dest="n_jobs")
parser.add_argument(
    "--no-bras",
    action="store_true",
    help="Disable silhouette batch (BRAS); avoids failures on tiny or degenerate batch×label data.",
)
args_ns, _ = parser.parse_known_args(shlex.split(args_str) if args_str.strip() else [])

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

if (adata.obs["label"].astype(str) == "Unknown").all():
    warnings.warn(
        "All cells have label 'Unknown'; bio-conservation metrics are not meaningful."
    )

ad_bm = adata.copy()
if "counts" not in ad_bm.layers:
    ad_bm.layers["counts"] = ad_bm.X.copy()
ad_bm.X = ad_bm.layers["counts"].copy()
sc.pp.normalize_total(ad_bm, target_sum=1e4)
sc.pp.log1p(ad_bm)

bm_kw = {}
if args_ns.no_bras:
    bm_kw["batch_correction_metrics"] = replace(BatchCorrection(), bras=False)

bm = Benchmarker(
    ad_bm,
    batch_key="batch",
    label_key="label",
    embedding_obsm_keys=["X_emb"],
    pre_integrated_embedding_obsm_key=None,
    n_jobs=args_ns.n_jobs,
    progress_bar=False,
    **bm_kw,
)
bm.prepare()
bm.benchmark()
results = bm.get_results(min_max_scale=False, clean_names=True)
results.to_csv(f"{prefix}_{integration}_metrics.tsv", sep="\t")


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


# MultiQC custom table: data is { sample (row): { column: value, ... } }.
# Fixed section id merges all integration runs (one row each); one embedding row per run.
metrics_row = results.iloc[0]
row = {str(c): _mqc_table_cell(metrics_row[c]) for c in results.columns}
mqc_data = {integration: row}
mqc_headers = {str(c): {"format": "{:.3f}"} for c in results.columns}

with open(f"{prefix}_{integration}_mqc.json", "w") as f_json:
    json.dump(
        {
            "id": "scib_metrics_benchmark",
            "plot_type": "table",
            "section_name": "scib-metrics",
            "description": "scib-metrics benchmark (one row per integration method).",
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
