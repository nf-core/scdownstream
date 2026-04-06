#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import json
import platform
import base64
import pickle
import argparse
import shlex

import numpy as np

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
method = "${method}"
args = "${args}"
parser = argparse.ArgumentParser()
parser.add_argument("--decimals", type=int, default=None)
params = parser.parse_args(shlex.split(args))

filter_col = "${filter_col ?: ''}"
filter_val = "${filter_val ?: ''}"

meta_id = "${meta.id}"
obs_key = "${obs_key}"

if filter_col and filter_val:
    adata = adata[adata.obs[filter_col] == filter_val].copy()

kwargs = {
    "groupby": obs_key,
    "method": method,
    "pts": True
}

# Check value counts for each group
value_counts = adata.obs[obs_key].value_counts()
# Filter out groups with less than 2 samples (scanpy requirement)
valid_groups = value_counts[value_counts >= 2].index.tolist()
invalid_groups = value_counts[value_counts < 2].index.tolist()

if len(invalid_groups) > 0:
    print(f"Warning: Excluding groups with < 2 samples: {', '.join(map(str, invalid_groups))}")

# Only proceed if we have at least 2 groups with >= 2 samples each
if len(valid_groups) >= 2:
    # Filter adata to only include valid groups
    adata = adata[adata.obs[obs_key].isin(valid_groups)].copy()

    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, **kwargs)

    rgg_dict = adata.uns["rank_genes_groups"]

    if params.decimals is not None:
        for key, arr in rgg_dict.items():
            if not hasattr(arr, 'dtype'):
                continue
            if np.issubdtype(arr.dtype, np.floating):
                rgg_dict[key] = arr.round(params.decimals)
            elif arr.dtype.names:
                rounded = np.empty_like(arr)
                for field in arr.dtype.names:
                    if np.issubdtype(arr.dtype[field], np.floating):
                        rounded[field] = np.round(arr[field], params.decimals)
                    else:
                        rounded[field] = arr[field]
                rgg_dict[key] = rounded

        # Re-sort deterministically: score desc, gene name asc to break ties
        groups = rgg_dict['names'].dtype.names
        sort_indices = {
            g: np.lexsort((
                np.array(rgg_dict['names'][g], dtype=str),
                -rgg_dict['scores'][g].astype(float)
            ))
            for g in groups
        }
        for key, arr in rgg_dict.items():
            if not hasattr(arr, 'dtype') or not arr.dtype.names:
                continue
            reordered = np.empty_like(arr)
            for g in arr.dtype.names:
                reordered[g] = arr[g][sort_indices[g]]
            rgg_dict[key] = reordered

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

        # Build section name with filter and obs_key information
        if filter_col and filter_val:
            section_name = f"Characteristic genes (grouped by: {obs_key}, filtered: {filter_col}={filter_val})"
            description = f"Characteristic genes, grouped by <code>{obs_key}</code>, filtered to <code>{filter_col}={filter_val}</code>."
        else:
            section_name = f"Characteristic genes (grouped by: {obs_key})"
            description = f"Characteristic genes, grouped by <code>{obs_key}</code>."

        custom_json = {
            "id": "${prefix}",
            "parent_id": "${meta.integration}",
            "parent_name": "${meta.integration}",
            "parent_description": "Results of the ${meta.integration} integration.",

            "section_name": section_name,
            "description": description,
            "plot_type": "image",
            "data": image_html,
        }

        json.dump(custom_json, f_json)
else:
    if len(valid_groups) == 0:
        print("Skipping rank_genes_groups computation: no groups have >= 2 samples.")
    elif len(valid_groups) == 1:
        print(f"Skipping rank_genes_groups computation: only one group has >= 2 samples (group: {valid_groups[0]}).")
    else:
        print("Skipping rank_genes_groups computation: less than 2 valid groups remaining after filtering.")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
