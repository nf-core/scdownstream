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
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"
mito_genes = "${mito_genes}"

log1p_total_counts_nmads = float("${log1p_total_counts_nmads}")
log1p_n_genes_by_counts_nmads = float("${log1p_n_genes_by_counts_nmads}")
pct_counts_in_top_20_genes_nmads = float("${pct_counts_in_top_20_genes_nmads}")
pct_counts_mt_nmads = float("${pct_counts_mt_nmads}")

mad_enabled = any(
    nmads > 0
    for nmads in [
        log1p_total_counts_nmads,
        log1p_n_genes_by_counts_nmads,
        pct_counts_in_top_20_genes_nmads,
        pct_counts_mt_nmads,
    ]
)

symbols = adata.var_names if symbol_col == "index" else adata.var[symbol_col]

if mito_genes:
    with open(mito_genes, "r") as f:
        mito_genes = {
            line.strip().lower()
            for line in f
            if line.strip() and not line.startswith("#")
        }
    adata.var["mt"] = symbols.str.lower().isin(mito_genes)
else:
    adata.var["mt"] = symbols.str.lower().str.startswith("mt-")

adata.var["ribo"] = adata.var_names.str.lower().str.match(r"^rp[sl]")
adata.var["hb"] = adata.var_names.str.lower().str.match(r"^hb[^p]")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=[20] if mad_enabled else None,
    log1p=mad_enabled,
    inplace=True,
)

if mad_enabled:
    mad_outlier = np.zeros(adata.n_obs, dtype=bool)

    if log1p_total_counts_nmads > 0:
        mad_outlier |= is_outlier(adata, "log1p_total_counts", log1p_total_counts_nmads)
    if log1p_n_genes_by_counts_nmads > 0:
        mad_outlier |= is_outlier(
            adata, "log1p_n_genes_by_counts", log1p_n_genes_by_counts_nmads
        )
    if pct_counts_in_top_20_genes_nmads > 0:
        mad_outlier |= is_outlier(
            adata, "pct_counts_in_top_20_genes", pct_counts_in_top_20_genes_nmads
        )
    if pct_counts_mt_nmads > 0:
        mad_outlier |= is_outlier(adata, "pct_counts_mt", pct_counts_mt_nmads)

    print(f"Total number of cells before MAD filtering: {adata.n_obs}")
    adata = adata[~mad_outlier].copy()
    print(f"Number of cells after MAD filtering: {adata.n_obs}")

adata = adata[adata.obs.pct_counts_mt < int("${max_mito_percentage}"), :].copy()
adata = adata[adata.obs.pct_counts_ribo >= int("${min_ribo_percentage}"), :].copy()
adata = adata[adata.obs.pct_counts_hb < int("${max_hb_percentage}"), :].copy()

sc.pp.filter_cells(adata, min_counts=int("${min_counts_cell}"))
sc.pp.filter_genes(adata, min_counts=int("${min_counts_gene}"))

sc.pp.filter_cells(adata, min_genes=int("${min_genes}"))
sc.pp.filter_genes(adata, min_cells=int("${min_cells}"))

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
