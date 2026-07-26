#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import base64
import json
import platform

import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import yaml

adata = sc.read_h5ad("${h5ad}")
symbol_col = "${symbol_col}"
mito_genes = "${mito_genes}"

if symbol_col in ("index", "none", ""):
    symbols = adata.var_names
else:
    symbols = adata.var[symbol_col]

if mito_genes:
    with open(mito_genes) as f:
        mito_genes = {line.strip().lower() for line in f if line.strip() and not line.startswith("#")}
    adata.var["mt"] = symbols.str.lower().isin(mito_genes)
else:
    adata.var["mt"] = symbols.str.lower().str.startswith("mt-")

has_mito = adata.var["mt"].any()

scatter_kwargs = {"x": "total_counts", "y": "n_genes_by_counts", "show": False}
if has_mito:
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    scatter_kwargs["color"] = "pct_counts_mt"
else:
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

sc.pl.scatter(adata, **scatter_kwargs)
path = "${prefix}_total_counts_vs_n_genes_by_counts.png"
plt.savefig(path, bbox_inches="tight")

# MultiQC

with open(path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
    image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    custom_json = {
        "id": "${prefix}",
        "parent_id": "${section_name}".replace(" ", "_"),
        "parent_name": "${section_name}",
        "parent_description": "${description}",
        "section_name": "${meta.id}",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "matplotlib": matplotlib.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
