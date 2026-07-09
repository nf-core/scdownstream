#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import base64
import json
import pickle
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
method = "${method}"
rank_key = "${rank_key}"

filter_col = "${filter_col ?: ''}"
filter_val = "${filter_val ?: ''}"

meta_id = "${meta.id}"
obs_key = "${obs_key}"

if filter_col and filter_val:
    adata = adata[adata.obs[filter_col] == filter_val].copy()

kwargs = {"groupby": obs_key, "method": method, "pts": True, "key_added": rank_key}

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
    sc.tl.filter_rank_genes_groups(
        adata,
        key=rank_key,
        min_in_group_fraction=0.2,
        max_out_group_fraction=0.2,
    )

    rgg_dict = adata.uns[rank_key]

    marker_df = sc.get.rank_genes_groups_df(adata, key=rank_key)
    marker_df.to_csv(f"{prefix}_markers.csv", index=False)

    pickle.dump(rgg_dict, open(f"{prefix}.pkl", "wb"))
    adata.write_h5ad(f"{prefix}.h5ad")

    # Plot
    sc.pl.rank_genes_groups(adata, key=rank_key, show=False)
    path = f"{prefix}.png"
    plt.savefig(path)

    sc.pl.rank_genes_groups_dotplot(adata, key=rank_key, show=False)
    dotplot_path = f"{prefix}_dotplot.png"
    plt.savefig(dotplot_path)

    # Build section name with filter and obs_key information
    if filter_col and filter_val:
        section_name = f"Characteristic genes (grouped by: {obs_key}, filtered: {filter_col}={filter_val})"
        description = f"Characteristic genes, grouped by <code>{obs_key}</code>, filtered to <code>{filter_col}={filter_val}</code>."
    else:
        section_name = f"Characteristic genes (grouped by: {obs_key})"
        description = f"Characteristic genes, grouped by <code>{obs_key}</code>."

    def write_mqc_plot(plot_path, plot_id, plot_label):
        with open(plot_path, "rb") as f_plot:
            image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        custom_json = {
            "id": plot_id,
            "parent_id": "${meta.integration}",
            "parent_name": "${meta.integration}",
            "parent_description": "Results of the ${meta.integration} integration.",
            "section_name": f"{section_name} ({plot_label})",
            "description": f"{description} {plot_label.capitalize()}.",
            "plot_type": "image",
            "data": image_html,
        }
        with open(f"{plot_id}_mqc.json", "w") as f_json:
            json.dump(custom_json, f_json)

    write_mqc_plot(path, "${prefix}", "rank plot")
    write_mqc_plot(dotplot_path, "${prefix}_dotplot", "dot plot")
else:
    if len(valid_groups) == 0:
        print("Skipping rank_genes_groups computation: no groups have >= 2 samples.")
    elif len(valid_groups) == 1:
        print(f"Skipping rank_genes_groups computation: only one group has >= 2 samples (group: {valid_groups[0]}).")
    else:
        print("Skipping rank_genes_groups computation: less than 2 valid groups remaining after filtering.")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
