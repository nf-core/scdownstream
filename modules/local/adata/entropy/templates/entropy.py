#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import json
import base64
import yaml
import argparse
import shlex

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
from scipy.stats import entropy
import scipy
import numpy as np
import matplotlib.pyplot as plt

group_col = "${group_col}"
entropy_col = "${entropy_col}"
prefix = "${prefix}"
adata = sc.read_h5ad("${h5ad}", backed='r')
args = "${args}"

parser = argparse.ArgumentParser()
parser.add_argument("--decimals", type=int, default=None)
params = parser.parse_args(shlex.split(args))

def entropy_of_group(group):
    counts = group.value_counts(normalize=True)
    return entropy(counts, base=2)

entropies = adata.obs.groupby(group_col)[entropy_col].apply(entropy_of_group)

n_unique = adata.obs[entropy_col].nunique()
colname = "${meta.id}:entropy"
adata.obs[colname] = adata.obs[group_col].map(entropies).astype(float) / np.log2(n_unique)

if params.decimals is not None:
    adata.obs[colname] = adata.obs[colname].round(params.decimals)

adata.obs[[colname]].to_pickle(f"{prefix}.pkl")
adata.write_h5ad(f"{prefix}.h5ad")

if "${plot_basis ? 'true' : 'false'}" == "true":
    sc.pl.scatter(adata, title="${meta.id} Entropy", color=colname, show=False, basis="${plot_basis}")
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

            "section_name": "${meta.id} Entropy",
            "plot_type": "image",
            "data": image_html,
        }

        json.dump(custom_json, f_json)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "scipy": scipy.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
