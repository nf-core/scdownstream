#!/usr/bin/env python3

import base64
import importlib.metadata
import json
import platform
from pathlib import Path

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
methods = "${methods.join(' ')}".split()
prediction_paths = r"""${predictions.join('\n')}""".splitlines()


def load_predictions(path):
    frame = pd.read_parquet(path, engine="pyarrow").set_index("__index__").rename_axis(None)
    if len(frame.columns) != 1:
        raise ValueError(f"{path} must contain exactly one Boolean column")
    if not frame.index.is_unique or set(frame.index) != set(adata.obs_names):
        raise ValueError(f"{path} must contain each input barcode exactly once")
    calls = frame.iloc[:, 0]
    if not pd.api.types.is_bool_dtype(calls) or calls.isna().any():
        raise ValueError(f"{path} must contain non-null Boolean predictions")
    return calls.reindex(adata.obs_names)


if not adata.obs_names.is_unique:
    raise ValueError("Input H5AD contains duplicate barcodes")

votes = pd.DataFrame(index=adata.obs_names)
annotations = pd.DataFrame(index=adata.obs_names)
for method, prediction_path in zip(methods, prediction_paths, strict=True):
    vote = load_predictions(prediction_path)
    votes[method] = vote
    annotations[f"{method}_predicted_doublet"] = vote
    annotations[f"{method}_class"] = vote.map({False: "singlet", True: "doublet"})

annotations = annotations.reindex(sorted(annotations.columns), axis=1)

if remove_doublets:
    adata = adata[~(votes.sum(axis=1) >= threshold), :].copy()

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")
annotations.reset_index(names="__index__").to_parquet(
    Path(f"{prefix}.obs.parquet").resolve(), index=False, engine="pyarrow", compression="snappy"
)

if len(votes.columns) > 1:
    contents = {column: votes.index[votes[column]].tolist() for column in votes}
    if any(contents.values()):
        plot_data = upsetplot.from_contents(contents)
        intersection_sizes = plot_data.groupby(level=list(range(len(contents)))).size()
        # The usual cutoff can hide every intersection in small samples.
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
    plt.close()
    with open(plot_path, "rb") as handle:
        image_string = base64.b64encode(handle.read()).decode("utf-8")
    custom_json = {
        "id": "${prefix}",
        "parent_id": "doublet_predictions",
        "parent_name": "Doublet predictions",
        "parent_description": "Upset plots of the various doublet prediction tools for each sample.",
        "section_name": "${meta.id}",
        "plot_type": "image",
        "data": f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>',
    }
    with open(f"{prefix}_mqc.json", "w") as handle:
        json.dump(custom_json, handle)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "pyarrow": importlib.metadata.version("pyarrow"),
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__,
    }
}
with open("versions.yml", "w") as handle:
    yaml.dump(versions, handle)
