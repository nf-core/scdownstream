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

symbol_col = "${symbol_col}"
prefix     = "${prefix}"


def read_gene_list(path):
    with open(path) as f:
        return [
            line.strip()
            for line in f
            if line.strip() and not line.startswith("#")
        ]


s_genes   = read_gene_list("${s_genes}")
g2m_genes = read_gene_list("${g2m_genes}")

adata = sc.read_h5ad("${h5ad}")

# If gene symbols are in a column rather than the index, temporarily set the
# index so that sc.tl.score_genes_cell_cycle can find them, then restore.
original_index = None
if symbol_col != "index":
    if symbol_col not in adata.var.columns:
        raise ValueError(
            f"symbol_col '{symbol_col}' not found in adata.var. "
            f"Available columns: {list(adata.var.columns)}"
        )
    original_index = adata.var_names.copy()
    adata.var_names = adata.var[symbol_col].astype(str)

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

if original_index is not None:
    adata.var_names = original_index

adata.obs[["S_score", "G2M_score", "phase"]].to_pickle(f"{prefix}.pkl")
adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
