#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"

import anndata as ad
from scimilarity.utils import pseudobulk_anndata
import scimilarity
import yaml

adata = ad.read_h5ad("${h5ad}")

counts_layer = "${counts_layer}"
groupby_labels = "${groupby_labels.join(',')}".split(',')
min_num_cells = int("${min_num_cells}")

assert counts_layer == "X" or counts_layer in adata.layers, f"Counts layer {counts_layer} not found in adata.layers"
if counts_layer != "counts":
    adata.layers["counts"] = adata.layers[counts_layer] if counts_layer != "X" else adata.X

adata_pseudobulk = pseudobulk_anndata(
    adata,
    groupby_labels,
    min_num_cells=min_num_cells
)

# Write the output
adata_pseudobulk.write_h5ad("${prefix}.h5ad")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scimilarity": scimilarity.__version__,
        "anndata": ad.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
