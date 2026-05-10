#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import platform
import base64
import json
import yaml

adata = sc.read_h5ad("${h5ad}")

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False)
path = "${prefix}_total_counts_vs_n_genes_by_counts.png"
plt.savefig(path)

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
