#!/usr/bin/env python3

import base64
import json
import platform
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy
import yaml
from scipy.stats import entropy


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        # An absolute path keeps zarr from parsing colons in Nextflow prefixes as URL schemes
        return ad.read_zarr(Path(path).absolute())
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(Path(f"{stem}.zarr").absolute())


group_col = "${group_col}"
entropy_col = "${entropy_col}"
prefix = "${prefix}"
adata = read_anndata("${anndata}", backed="r")


def entropy_of_group(group):
    counts = group.value_counts(normalize=True)
    return entropy(counts, base=2)


entropies = adata.obs.groupby(group_col)[entropy_col].apply(entropy_of_group)

n_unique = adata.obs[entropy_col].nunique()
colname = "${meta.id}:entropy"
adata.obs[colname] = adata.obs[group_col].map(entropies).astype(float) / np.log2(n_unique)

adata.obs[[colname]].to_pickle(f"{prefix}.pkl")
write_anndata(adata, prefix, "${anndata_format}")

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
