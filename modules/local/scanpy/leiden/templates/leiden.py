#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import json
import base64
import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import matplotlib.pyplot as plt

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}", backed='r')
prefix = "${prefix}"
key_added = "${key_added}"

kwargs = {
    "resolution": float("${resolution}"),
    "key_added": key_added
}

sc.tl.leiden(adata, **kwargs)

adata.obs[[key_added]].to_pickle(f"{prefix}.pkl")
adata.write_h5ad(f"{prefix}.h5ad")

if "${plot_umap}" == "true":
    # Plot
    sc.pl.umap(adata, title="${meta.id} Leiden", color=key_added, show=False)
    path = f"{prefix}.png"
    plt.savefig(path, bbox_inches='tight')

    # MultiQC
    with open(path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

        custom_json = {
            "id": "${prefix}",
            "parent_id": "${meta.integration}",
            "parent_name": "${meta.integration}",
            "parent_description": "Results of the ${meta.integration} integration.",

            "section_name": "${meta.id} Leiden",
            "plot_type": "image",
            "data": image_html,
        }

        json.dump(custom_json, f_json)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
