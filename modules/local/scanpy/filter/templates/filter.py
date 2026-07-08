#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import numpy as np
import scanpy as sc
import yaml
from scipy.stats import median_abs_deviation
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


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


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"
mito_genes = "${mito_genes}"

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

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
