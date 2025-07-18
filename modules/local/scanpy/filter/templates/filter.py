#!/usr/bin/env python3

import os
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


def is_outlier(adata, metric: str, nmads: int):
    """Identify outliers using MAD (Median Absolute Deviation) method."""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True
)

if "${dynamic_filtering}" == "true":
    # Dynamic filtering using MAD strategy
    print("Applying dynamic filtering using MAD strategy...")

    # Filter cells based on general QC metrics (5 MADs)
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )

    # Filter mitochondrial outliers (3 MADs) with additional max threshold
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > int("${max_mito_percentage}")
    )

    # Apply filtering
    print(f"Total number of cells before filtering: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after dynamic filtering: {adata.n_obs}")

else:
    # Static filtering using user-defined thresholds
    print("Applying static filtering using user-defined thresholds...")

    # Apply mitochondrial filtering
    adata = adata[adata.obs.pct_counts_mt < int("${max_mito_percentage}"), :].copy()

    # Apply count and gene filtering
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
