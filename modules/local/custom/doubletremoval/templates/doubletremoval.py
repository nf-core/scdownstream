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

h5ad_path = "${h5ad}"
adata = ad.read_h5ad(h5ad_path)
threshold = int("${threshold}")
remove_doublets = "${removal}" == "true"
prefix = "${prefix}"
# Detector modules name their call column "<sample id>_<method>". Dropping the sample id
# gives column names that are identical in every dataset, so they survive merging.
sample_prefix = "${meta.id}_"


def load(path: str) -> pd.DataFrame:
    if path.endswith(".pkl"):
        return pd.read_pickle(path)
    if path.endswith(".csv"):
        return pd.read_csv(path, index_col=0)
    raise ValueError(f"Unsupported doublet prediction table '{path}': expected a .csv or .pkl file.")


def examples(values) -> str:
    values = sorted(str(value) for value in values)
    return ", ".join(values[:5]) + (", ..." if len(values) > 5 else "")


def is_call(column: pd.Series) -> bool:
    """Boolean columns are doublet calls (votes); everything else is an annotation."""
    if pd.api.types.is_bool_dtype(column):
        return True
    values = column.dropna()
    return column.dtype == object and len(values) > 0 and values.map(lambda value: isinstance(value, bool)).all()


if not adata.obs_names.is_unique:
    duplicated = adata.obs_names[adata.obs_names.duplicated()].unique()
    raise ValueError(
        f"AnnData '{h5ad_path}' has {len(duplicated)} duplicated cell barcodes ({examples(duplicated)}); "
        "doublet predictions cannot be joined unambiguously."
    )

calls = {}
annotations = {}

for path in sorted("${predictions}".split()):
    table = load(path)
    table.index = table.index.astype(str)
    table.columns = [str(column).removeprefix(sample_prefix) for column in table.columns]

    # Join on the exact cell barcode. Cells without a prediction must not silently
    # become singlets, and predictions for unknown cells point to a mismatched input.
    if not table.index.is_unique:
        duplicated = table.index[table.index.duplicated()].unique()
        raise ValueError(
            f"Doublet prediction table '{path}' has {len(duplicated)} duplicated cell barcodes ({examples(duplicated)})."
        )
    missing = adata.obs_names.difference(table.index)
    unexpected = table.index.difference(adata.obs_names)
    if len(missing) or len(unexpected):
        raise ValueError(
            f"Doublet prediction table '{path}' does not match the cells in '{h5ad_path}': "
            f"{len(missing)} cells have no prediction ({examples(missing)}) and "
            f"{len(unexpected)} predictions have no matching cell ({examples(unexpected)})."
        )
    table = table.loc[adata.obs_names]

    call_columns = [column for column in table.columns if is_call(table[column])]
    if not call_columns:
        raise ValueError(f"Doublet prediction table '{path}' has no boolean doublet call column.")

    for column in table.columns:
        if column in calls or column in annotations:
            raise ValueError(f"Doublet prediction column '{column}' is provided by more than one prediction table.")
        if column in call_columns:
            if table[column].isna().any():
                raise ValueError(
                    f"Doublet call column '{column}' in '{path}' has {int(table[column].isna().sum())} missing values."
                )
            calls[column] = table[column].astype(bool)
        else:
            # Annotations stay in the namespace of a call column from the same table,
            # so that different methods can never overwrite each other.
            if not any(column.startswith(f"{call}_") for call in call_columns):
                raise ValueError(
                    f"Annotation column '{column}' in '{path}' must start with one of: "
                    + ", ".join(f"'{call}_'" for call in call_columns)
                )
            values = table[column]
            if not pd.api.types.is_numeric_dtype(values):
                values = values.astype("category")
            annotations[column] = values

# Only the boolean calls take part in the vote; scores and classes are annotations.
# Stable column order so nf-test snapshots do not depend on process finish order.
calls = pd.DataFrame(calls, index=adata.obs_names)
calls = calls.reindex(sorted(calls.columns), axis=1)
annotations = pd.DataFrame(annotations, index=adata.obs_names)

doublet_columns = pd.concat([calls, annotations], axis=1)
for column in sorted(doublet_columns.columns):
    if column in adata.obs.columns:
        print(f"Overwriting existing obs column '{column}' with doublet detection results.")
    adata.obs[column] = doublet_columns[column]

if remove_doublets:
    mask = (calls.sum(axis=1) >= threshold).to_numpy()
    adata = adata[~mask, :].copy()

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

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

if not len(calls.columns) > 1:
    score_columns = [column for column in calls.columns if "score" in column.lower()]
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

contents = {column: calls[column][calls[column]].index.tolist() for column in calls.columns}

if any(contents.values()):
    plot_data = upsetplot.from_contents(contents)
    intersection_sizes = plot_data.groupby(level=list(range(len(contents)))).size()
    # Keep the normal display cutoff, unless it would hide every intersection.
    minimum_size = 10 if intersection_sizes.max() >= 10 else 1
    upsetplot.plot(
        plot_data,
        sort_by="cardinality",
        show_counts=True,
        subset_size="count",
        min_subset_size=minimum_size,
    )
else:
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.text(0.5, 0.5, "No doublets predicted", ha="center", va="center", transform=ax.transAxes)
    ax.set_axis_off()
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
