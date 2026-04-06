#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import yaml
import scanpy as sc
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"
mito_genes = "${mito_genes}"

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

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
adata = adata[adata.obs.pct_counts_mt < int("${max_mito_percentage}"), :].copy()

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
