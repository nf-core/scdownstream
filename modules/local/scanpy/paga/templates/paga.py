#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import json
import base64
import pickle
import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}", backed='r')
prefix = "${prefix}"
obs_key = "${obs_key}"

if adata.obs[obs_key].value_counts().size > 1:
    sc.tl.paga(adata, groups=obs_key)

    paga_dict = adata.uns["paga"]

    # Save PAGA data
    pickle.dump(paga_dict, open(f"{prefix}.pkl", "wb"))

    np.save(f"{prefix}_connectivities.npy", adata.obsp["connectivities"])
    adata.write_h5ad(f"{prefix}.h5ad")

    # Plot
    sc.pl.paga(adata, title="${meta.id} PAGA", show=False)
    path = f"{prefix}.png"
    plt.savefig(path)

    # MultiQC
    with open(path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

        custom_json = {
            "id": "${prefix}",
            "parent_id": "${meta.integration}",
            "parent_name": "${meta.integration}",
            "parent_description": "Results of the ${meta.integration} integration.",

            "section_name": "${meta.id} PAGA",
            "plot_type": "image",
            "data": image_html,
        }

        json.dump(custom_json, f_json)
else:
    print(f"Skipping PAGA computation for {obs_key} as it has less than 2 unique values.")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
