#!/usr/bin/env python3

import base64
import json
import platform

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot
import yaml


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        return ad.read_zarr(path)
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(f"{stem}.zarr")


adata = read_anndata("${anndata}")
threshold = int("${threshold}")
remove_doublets = "${removal}" == "true"
prefix = "${prefix}"


def load(path: str) -> pd.DataFrame:
    if path.endswith(".pkl"):
        return pd.read_pickle(path)
    if path.endswith(".csv"):
        return pd.read_csv(path, index_col=0)


predictions = pd.concat([load(f) for f in "${predictions}".split()], axis=1)
predictions = predictions.reindex(adata.obs_names)
# Stable column order so nf-test snapshots do not depend on process finish order.
predictions = predictions.reindex(sorted(predictions.columns), axis=1)
for column in predictions.columns:
    adata.obs[column] = predictions[column]

if remove_doublets:
    mask = predictions.sum(axis=1) >= threshold
    adata = adata[~mask, :]

write_anndata(adata, prefix, "${anndata_format}")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)

if not len(predictions.columns) > 1:
    score_columns = [column for column in predictions.columns if "score" in column.lower()]
    if score_columns:
        score_col = score_columns[0]
        fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        adata.obs[score_col].hist(ax=ax, bins=50)
        ax.set_xlabel(score_col)
        ax.set_ylabel("Cells")
        plot_path = f"{prefix}_doublet_scores.png"
        plt.savefig(plot_path)

        with open(plot_path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
            image_string = base64.b64encode(f_plot.read()).decode("utf-8")
            image_html = (
                f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
            )
            custom_json = {
                "id": "${prefix}",
                "parent_id": "doublet_predictions",
                "parent_name": "Doublet predictions",
                "parent_description": "Doublet score distributions and tool overlap per sample.",
                "section_name": "${meta.id}",
                "plot_type": "image",
                "data": image_html,
            }
            json.dump(custom_json, f_json)
    exit(0)

# Plot

contents = {column: predictions[column][predictions[column]].index.tolist() for column in predictions.columns}

plot_data = upsetplot.from_contents(contents)

upsetplot.plot(plot_data, sort_by="cardinality", show_counts=True, subset_size="count", min_subset_size=10)
plot_path = f"{prefix}_predictions_mqc.png"
plt.savefig(plot_path)

# MultiQC

with open(plot_path, "rb") as f_plot, open("${prefix}_mqc.json", "w") as f_json:
    image_string = base64.b64encode(f_plot.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    custom_json = {
        "id": "${prefix}",
        "parent_id": "doublet_predictions",
        "parent_name": "Doublet predictions",
        "parent_description": "Upset plots of the various doublet prediction tools for each sample.",
        "section_name": "${meta.id}",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)
