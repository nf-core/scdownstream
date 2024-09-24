#!/usr/bin/env python3

import json
import platform
import base64
import pickle

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

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


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
use_gpu = "${task.ext.use_gpu}" == "true"

kwargs = {
    "groupby": "${obs_key}",
    "pts": True
}


if use_gpu:
    import rapids_singlecell as rsc
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    import cupy as cp
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

    rsc.get.anndata_to_GPU(adata)
    rsc.pp.log1p(adata)
    rsc.tl.rank_genes_groups_logreg(adata, **kwargs)
    rsc.get.anndata_to_CPU(adata)
else:
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, **kwargs)

rgg_dict = adata.uns["rank_genes_groups"]

pickle.dump(rgg_dict, open(f"{prefix}.pkl", "wb"))
adata.write_h5ad(f"{prefix}.h5ad")

# Plot
sc.pl.rank_genes_groups(adata, show=False)
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

        "section_name": "${meta.id} characteristic genes",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
