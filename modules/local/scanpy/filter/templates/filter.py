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


def mad_bounds(values, nmads):
    median = np.median(values)
    mad = median_abs_deviation(values)
    return median - nmads * mad, median + nmads * mad


def nmads_for_metric(metric, thresholds):
    mad_nmads = {
        "log1p_total_counts": thresholds["log1p_total_counts_nmads"],
        "log1p_n_genes_by_counts": thresholds["log1p_n_genes_by_counts_nmads"],
        "pct_counts_in_top_20_genes": thresholds["pct_counts_in_top_20_genes_nmads"],
        "pct_counts_mt": thresholds["pct_counts_mt_nmads"],
    }
    return mad_nmads.get(metric)


LOWER_BOUND_COLOR = "#1f77b4"
UPPER_BOUND_COLOR = "#d62728"
KEPT_COLOR = "#4daf4a"
FILTERED_COLOR = "#e41a1c"
MAD_LINESTYLE = "--"
ABSOLUTE_LINESTYLE = ":"


def threshold_bounds(metric, thresholds):
    """Return all configured bounds to plot (MAD and absolute when both are set)."""
    values = thresholds["values"][metric]
    bounds = []
    nmads = nmads_for_metric(metric, thresholds)

    if nmads is not None and nmads > 0:
        mad_lower, mad_upper = mad_bounds(values, nmads)
        bounds.append({"value": mad_lower, "side": "lower", "source": "MAD"})
        bounds.append({"value": mad_upper, "side": "upper", "source": "MAD"})

    manual_lower = None
    manual_upper = None

    if metric == "pct_counts_mt":
        max_mito = thresholds["max_mito_percentage"]
        if max_mito is not None and max_mito < 100:
            manual_upper = float(max_mito)
    elif metric == "pct_counts_ribo":
        min_ribo = thresholds["min_ribo_percentage"]
        if min_ribo is not None and min_ribo > 0:
            manual_lower = float(min_ribo)
    elif metric == "pct_counts_hb":
        max_hb = thresholds["max_hb_percentage"]
        if max_hb is not None and max_hb < 100:
            manual_upper = float(max_hb)
    elif metric == "total_counts":
        min_counts_cell = thresholds["min_counts_cell"]
        if min_counts_cell is not None and min_counts_cell > 0:
            manual_lower = float(min_counts_cell)
    elif metric == "n_genes_by_counts":
        min_genes = thresholds["min_genes"]
        if min_genes is not None and min_genes > 0:
            manual_lower = float(min_genes)

    if manual_lower is not None:
        bounds.append({"value": manual_lower, "side": "lower", "source": "absolute"})
    if manual_upper is not None:
        bounds.append({"value": manual_upper, "side": "upper", "source": "absolute"})

    return bounds


def bound_legend_label(side, source):
    bound_side = "Lower" if side == "lower" else "Upper"
    bound_type = "MAD" if source == "MAD" else "absolute"
    return f"{bound_side} bound ({bound_type})"


def cell_keep_mask(adata, mad_filters):
    """Whether each cell is retained after all cell-removing filter steps."""
    keep = np.ones(adata.n_obs, dtype=bool)

    mad_outlier = np.zeros(adata.n_obs, dtype=bool)
    for metric, nmads in mad_filters:
        if nmads is not None and nmads > 0:
            mad_outlier |= is_outlier(adata, metric, nmads)
    keep &= ~mad_outlier

    if max_mito_percentage is not None:
        keep &= adata.obs.pct_counts_mt < max_mito_percentage
    if min_ribo_percentage is not None:
        keep &= adata.obs.pct_counts_ribo >= min_ribo_percentage
    if max_hb_percentage is not None:
        keep &= adata.obs.pct_counts_hb < max_hb_percentage
    if min_counts_cell is not None:
        keep &= adata.obs.total_counts >= min_counts_cell
    if min_genes is not None:
        keep &= adata.obs.n_genes_by_counts >= min_genes

    return keep


def plot_metric_histogram(ax, metric, values, cell_keep, thresholds):
    """Draw a QC metric histogram on the given axes."""
    bins = np.histogram_bin_edges(values, bins=50)
    values_kept = values[cell_keep]
    values_filtered = values[~cell_keep]

    hist_series = []
    hist_colors = []
    hist_labels = []
    if values_kept.size > 0:
        hist_series.append(values_kept)
        hist_colors.append(KEPT_COLOR)
        hist_labels.append("Kept")
    if values_filtered.size > 0:
        hist_series.append(values_filtered)
        hist_colors.append(FILTERED_COLOR)
        hist_labels.append("Filtered out")

    ax.hist(
        hist_series,
        bins=bins,
        stacked=len(hist_series) > 1,
        color=hist_colors,
        edgecolor="white",
        linewidth=0.5,
        label=hist_labels,
    )

    median = np.median(values)
    ax.axvline(median, color="black", linestyle="-.", linewidth=1.2, label="Median")

    for bound in threshold_bounds(metric, thresholds):
        side = bound["side"]
        source = bound["source"]
        color = LOWER_BOUND_COLOR if side == "lower" else UPPER_BOUND_COLOR
        linestyle = MAD_LINESTYLE if source == "MAD" else ABSOLUTE_LINESTYLE
        ax.axvline(
            bound["value"],
            color=color,
            linestyle=linestyle,
            linewidth=1.2,
            label=bound_legend_label(side, source),
        )

    ax.set_xlabel(metric, fontsize=8)
    ax.set_ylabel("Cells", fontsize=8)
    ax.set_title(metric, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.legend(loc="upper right", fontsize=6)


def plot_qc_histogram_panel(prefix, section_name, description, cell_keep, thresholds):
    """Plot all QC metric histograms on a 3x3 grid (bottom-right panel empty)."""
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    axes = axes.flatten()

    for index, metric in enumerate(PLOT_METRICS):
        plot_metric_histogram(
            axes[index],
            metric,
            thresholds["values"][metric],
            cell_keep,
            thresholds,
        )

    axes[8].set_visible(False)
    fig.tight_layout()

    path = f"{prefix}_qc_histograms.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)

    with open(path, "rb") as f_plot, open(f"{prefix}_qc_histograms_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

        custom_json = {
            "id": f"{prefix}_qc_histograms",
            "parent_id": section_name.replace(" ", "_"),
            "parent_name": section_name,
            "parent_description": description,
            "section_name": "${meta.id}",
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
    cell_keep = cell_keep_mask(adata, mad_filters)
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
    plot_qc_histogram_panel(
        prefix,
        section_name,
        description,
        cell_keep,
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
