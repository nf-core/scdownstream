#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import base64
import json
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import yaml
from scipy.stats import median_abs_deviation
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

PLOT_METRICS = [
    "log1p_total_counts",
    "log1p_n_genes_by_counts",
    "pct_counts_in_top_20_genes",
    "pct_counts_mt",
    "pct_counts_ribo",
    "pct_counts_hb",
    "total_counts",
    "n_genes_by_counts",
]


def is_outlier(adata, metric: str, nmads: float):
    """Identify outliers using MAD (Median Absolute Deviation) method."""
    metric_values = adata.obs[metric]
    outlier = (metric_values < np.median(metric_values) - nmads * median_abs_deviation(metric_values)) | (
        np.median(metric_values) + nmads * median_abs_deviation(metric_values) < metric_values
    )
    return outlier


def parse_optional_number(value):
    if value is None:
        return None
    text = str(value).strip()
    if text in ("", "null", "None"):
        return None
    return float(text)


def parse_optional_int(value):
    parsed = parse_optional_number(value)
    return None if parsed is None else int(parsed)


def mad_thresholds(values, nmads):
    median = np.median(values)
    mad = median_abs_deviation(values)
    return [median - nmads * mad, median + nmads * mad]


def threshold_lines(metric, thresholds):
    """Return dashed threshold lines for a metric."""
    lines = []

    mad_nmads = {
        "log1p_total_counts": thresholds["log1p_total_counts_nmads"],
        "log1p_n_genes_by_counts": thresholds["log1p_n_genes_by_counts_nmads"],
        "pct_counts_in_top_20_genes": thresholds["pct_counts_in_top_20_genes_nmads"],
        "pct_counts_mt": thresholds["pct_counts_mt_nmads"],
    }
    if metric in mad_nmads:
        nmads = mad_nmads[metric]
        if nmads is not None and nmads > 0:
            lines.extend(mad_thresholds(thresholds["values"][metric], nmads))

    if metric == "pct_counts_mt":
        max_mito = thresholds["max_mito_percentage"]
        if max_mito is not None and max_mito < 100:
            lines.append(max_mito)
    elif metric == "pct_counts_ribo":
        min_ribo = thresholds["min_ribo_percentage"]
        if min_ribo is not None and min_ribo > 0:
            lines.append(min_ribo)
    elif metric == "pct_counts_hb":
        max_hb = thresholds["max_hb_percentage"]
        if max_hb is not None and max_hb < 100:
            lines.append(max_hb)
    elif metric == "total_counts":
        min_counts_cell = thresholds["min_counts_cell"]
        if min_counts_cell is not None and min_counts_cell > 0:
            lines.append(min_counts_cell)
    elif metric == "n_genes_by_counts":
        min_genes = thresholds["min_genes"]
        if min_genes is not None and min_genes > 0:
            lines.append(min_genes)

    return lines


def plot_qc_histogram(metric, values, prefix, section_name, description, thresholds):
    """Plot a QC metric histogram with median and threshold lines."""
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values, bins=50, color="steelblue", edgecolor="white", linewidth=0.5)

    median = np.median(values)
    ax.axvline(median, color="black", linestyle="-.", linewidth=1.5, label="Median")

    for threshold in threshold_lines(metric, thresholds):
        ax.axvline(threshold, color="crimson", linestyle="--", linewidth=1.5)

    ax.set_xlabel(metric)
    ax.set_ylabel("Cells")
    ax.set_title(metric)
    ax.legend(loc="upper right")

    path = f"{prefix}_{metric}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)

    with open(path, "rb") as f_plot, open(f"{prefix}_{metric}_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

        custom_json = {
            "id": f"{prefix}_{metric}",
            "parent_id": section_name.replace(" ", "_"),
            "parent_name": section_name,
            "parent_description": description,
            "section_name": "${meta.id} " + metric,
            "plot_type": "image",
            "data": image_html,
        }

        json.dump(custom_json, f_json)


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"
mito_genes = "${mito_genes}"
plot = "${plot}" == "true"
section_name = "${section_name}"
description = "${description}"

min_genes = parse_optional_int("${min_genes}")
min_cells = parse_optional_int("${min_cells}")
min_counts_gene = parse_optional_int("${min_counts_gene}")
min_counts_cell = parse_optional_int("${min_counts_cell}")
max_mito_percentage = parse_optional_int("${max_mito_percentage}")
min_ribo_percentage = parse_optional_int("${min_ribo_percentage}")
max_hb_percentage = parse_optional_int("${max_hb_percentage}")
log1p_total_counts_nmads = parse_optional_number("${log1p_total_counts_nmads}")
log1p_n_genes_by_counts_nmads = parse_optional_number("${log1p_n_genes_by_counts_nmads}")
pct_counts_in_top_20_genes_nmads = parse_optional_number("${pct_counts_in_top_20_genes_nmads}")
pct_counts_mt_nmads = parse_optional_number("${pct_counts_mt_nmads}")

mad_filters = [
    ("log1p_total_counts", log1p_total_counts_nmads),
    ("log1p_n_genes_by_counts", log1p_n_genes_by_counts_nmads),
    ("pct_counts_in_top_20_genes", pct_counts_in_top_20_genes_nmads),
    ("pct_counts_mt", pct_counts_mt_nmads),
]
mad_enabled = any(nmads is not None and nmads > 0 for _, nmads in mad_filters)

symbols = adata.var_names if symbol_col == "index" else adata.var[symbol_col]

if mito_genes:
    with open(mito_genes) as f:
        mito_genes = {line.strip().lower() for line in f if line.strip() and not line.startswith("#")}
    adata.var["mt"] = symbols.str.lower().isin(mito_genes)
else:
    adata.var["mt"] = symbols.str.lower().str.startswith("mt-")

adata.var["ribo"] = adata.var_names.str.lower().str.match(r"^rp[sl]")
adata.var["hb"] = adata.var_names.str.lower().str.match(r"^hb[^p]")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=[20],
    log1p=True,
    inplace=True,
)

if plot:
    threshold_context = {
        "values": {metric: adata.obs[metric].to_numpy() for metric in PLOT_METRICS},
        "log1p_total_counts_nmads": log1p_total_counts_nmads,
        "log1p_n_genes_by_counts_nmads": log1p_n_genes_by_counts_nmads,
        "pct_counts_in_top_20_genes_nmads": pct_counts_in_top_20_genes_nmads,
        "pct_counts_mt_nmads": pct_counts_mt_nmads,
        "max_mito_percentage": max_mito_percentage,
        "min_ribo_percentage": min_ribo_percentage,
        "max_hb_percentage": max_hb_percentage,
        "min_counts_cell": min_counts_cell,
        "min_genes": min_genes,
    }
    for metric in PLOT_METRICS:
        plot_qc_histogram(
            metric,
            threshold_context["values"][metric],
            prefix,
            section_name,
            description,
            threshold_context,
        )

if mad_enabled:
    mad_outlier = np.zeros(adata.n_obs, dtype=bool)

    for metric, nmads in mad_filters:
        if nmads is not None and nmads > 0:
            mad_outlier |= is_outlier(adata, metric, nmads)

    print(f"Total number of cells before MAD filtering: {adata.n_obs}")
    adata = adata[~mad_outlier].copy()
    print(f"Number of cells after MAD filtering: {adata.n_obs}")

if max_mito_percentage is not None:
    adata = adata[adata.obs.pct_counts_mt < max_mito_percentage, :].copy()
if min_ribo_percentage is not None:
    adata = adata[adata.obs.pct_counts_ribo >= min_ribo_percentage, :].copy()
if max_hb_percentage is not None:
    adata = adata[adata.obs.pct_counts_hb < max_hb_percentage, :].copy()

if min_counts_cell is not None:
    sc.pp.filter_cells(adata, min_counts=min_counts_cell)
if min_counts_gene is not None:
    sc.pp.filter_genes(adata, min_counts=min_counts_gene)
if min_genes is not None:
    sc.pp.filter_cells(adata, min_genes=min_genes)
if min_cells is not None:
    sc.pp.filter_genes(adata, min_cells=min_cells)

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}
if plot:
    versions["${task.process}"]["matplotlib"] = matplotlib.__version__

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
