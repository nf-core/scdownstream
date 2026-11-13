#!/usr/bin/env python3

import base64
import json
import platform

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy
import yaml
from scipy.stats import entropy

group_col = "${group_col}"
entropy_col = "${entropy_col}"
prefix = "${prefix}"
adata = sc.read_h5ad("${h5ad}", backed="r")

if group_col not in adata.obs.columns:
    raise ValueError(f"Grouping column '{group_col}' not found in adata.obs")
if entropy_col not in adata.obs.columns:
    raise ValueError(f"Entropy column '{entropy_col}' not found in adata.obs")

n_observed_labels = adata.obs[entropy_col].nunique(dropna=True)

if n_observed_labels == 0:
    raise ValueError(f"No observed non-missing labels in entropy column '{entropy_col}'")


def entropy_of_group(group):
    counts = group.value_counts(normalize=True)
    return entropy(counts, base=2)


colname = "${meta.id}:entropy"

if n_observed_labels == 1:
    adata.obs[colname] = 0.0
else:
    entropies = adata.obs.groupby(group_col)[entropy_col].apply(entropy_of_group)
    adata.obs[colname] = adata.obs[group_col].map(entropies).astype(float) / np.log2(n_observed_labels)

adata.obs[[colname]].to_pickle(f"{prefix}.pkl")
adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

if "${plot_basis ? 'true' : 'false'}" == "true":
    sc.pl.scatter(adata, title="${meta.id} Entropy", color=colname, show=False, basis="${plot_basis}")
    path = f"{prefix}.png"
    plt.savefig(path, bbox_inches="tight")

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
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "scipy": scipy.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
