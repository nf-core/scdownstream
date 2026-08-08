#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import base64
import json
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot
import yaml

adata = ad.read_h5ad("${h5ad}")
threshold = int("${threshold}")
remove_doublets = "${removal}" == "true"
prefix = "${prefix}"


def _dataframe_for_h5ad(df: pd.DataFrame) -> pd.DataFrame:
    """Rewrite string indexes/columns so nft-anndata can read the written H5AD.

    Newer anndata writes plain string indexes as nullable-string-array groups.
    nft-anndata treats every index group as categorical and NPEs without categories.
    Categorical indexes encode as categories/codes, which nft-anndata supports.
    """
    df = df.copy()
    df.index = pd.CategoricalIndex(df.index.astype(str).to_list(), name=df.index.name)
    for col in df.columns:
        if isinstance(df[col].dtype, pd.StringDtype) or df[col].dtype == object:
            df[col] = pd.Series(df[col].astype(str).to_list(), index=df.index, name=col, dtype=object)
    return df


def _prepare_adata_for_h5ad(adata_obj):
    """Avoid nullable-string-array encodings that nft-anndata misreads as categoricals."""
    adata_obj.obs = _dataframe_for_h5ad(adata_obj.obs)
    adata_obj.var = _dataframe_for_h5ad(adata_obj.var)
    for uns_key, uns_value in list(adata_obj.uns.items()):
        if isinstance(uns_value, pd.DataFrame):
            adata_obj.uns[uns_key] = _dataframe_for_h5ad(uns_value)
        elif isinstance(uns_value, dict):
            for key, value in list(uns_value.items()):
                if isinstance(value, pd.DataFrame):
                    uns_value[key] = _dataframe_for_h5ad(value)


def load(path: str) -> pd.DataFrame:
    if path.endswith(".parquet"):
        return pd.read_parquet(path)
    if path.endswith(".pkl"):
        return pd.read_pickle(path)
    if path.endswith(".csv"):
        return pd.read_csv(path, index_col=0)
    raise ValueError(f"Unsupported prediction file extension: {path}")


predictions = pd.concat([load(f) for f in "${predictions}".split()], axis=1)
predictions = predictions.reindex(adata.obs_names).fillna(False).astype(bool)
for column in predictions.columns:
    adata.obs[column] = predictions[column]

if remove_doublets:
    mask = predictions.sum(axis=1) >= threshold
    adata = adata[~mask, :]

_prepare_adata_for_h5ad(adata)
adata.write_h5ad(f"{prefix}.h5ad")

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

contents = {column: predictions.index[predictions[column]].tolist() for column in predictions.columns}

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
