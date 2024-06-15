#!/usr/bin/env python3

import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
import platform
import os
import scipy

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

adatas = {os.path.basename(f).split(".")[0]: sc.read_h5ad(f) for f in "${h5ads}".split()}
genes = [adata.var_names for adata in adatas.values()]

adata_outer = ad.concat(adatas, join="outer", index_unique="-")
adata_outer.X = csr_matrix(adata_outer.X)
adata_outer.layers["counts"] = adata_outer.X

# Make sure there are no cells and genes without any counts
sc.pp.filter_cells(adata_outer, min_counts=1)
sc.pp.filter_genes(adata_outer, min_cells=1)

genes_intersection = list(set(adata_outer.var_names).intersection(*genes))
sorted(genes_intersection)
adata_inner = adata_outer[:, genes_intersection]

adata_inner.write("${prefix}_inner.h5ad")
adata_outer.write("${prefix}_outer.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "scipy": scipy.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
