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

adata = ad.read_h5ad("${h5ad}")
threshold = int("${threshold}")
remove_doublets = "${removal}" == "true"
prefix = "${prefix}"
methods = json.loads("""${methods_json}""")
prediction_paths = json.loads("""${predictions_json}""")


def load(path, method):
    source = f"{method} table '{path}'"
    if path.endswith(".csv"):
        table = pd.read_csv(path, index_col=0, converters={0: str})
    elif path.endswith(".pkl"):
        table = pd.read_pickle(path)
    else:
        raise ValueError(f"{source} must be CSV or pickle")
    if len(table.columns) != 1:
        raise ValueError(f"{method} prediction '{path}' must contain exactly one Boolean column")
    if not table.index.is_unique:
        duplicates = table.index[table.index.duplicated()].tolist()[:5]
        raise ValueError(f"{source} contains duplicate indices: {duplicates}")
    missing = adata.obs_names.difference(table.index).tolist()
    extra = table.index.difference(adata.obs_names).tolist()
    if missing or extra:
        raise ValueError(f"{source} barcode coverage mismatch: missing={missing[:5]}, extra={extra[:5]}")
    return table.reindex(adata.obs_names)


def boolean_votes(values, source):
    # Explicit tokens avoid treating strings such as 'False' as positive votes.
    tokens = {"true": True, "false": False, "1": True, "0": False}
    normalised = values.map(lambda value: tokens.get(str(value).lower()))
    invalid = normalised.isna()
    if invalid.any():
        cells = values.index[invalid].tolist()[:5]
        raise ValueError(f"{source} contains null or invalid Boolean votes at {cells}")
    return normalised.astype(bool)


if not adata.obs_names.is_unique:
    duplicates = adata.obs_names[adata.obs_names.duplicated()].tolist()[:5]
    raise ValueError(f"Input H5AD contains duplicate barcodes: {duplicates}")
if not adata.obs.columns.is_unique:
    raise ValueError("Input H5AD contains duplicate observation columns")
supported = {"doubletdetection", "scdblfinder", "scrublet", "solo"}
if (
    not isinstance(methods, list)
    or not methods
    or any(not isinstance(method, str) or method not in supported for method in methods)
):
    raise ValueError(f"Unknown or empty doublet methods: {methods}")
if len(methods) != len(set(methods)):
    raise ValueError(f"Duplicate doublet methods: {methods}")
if methods != sorted(methods):
    raise ValueError(f"Doublet methods must be sorted: {methods}")
if len(prediction_paths) != len(methods):
    raise ValueError("Table counts must match methods: one prediction per method")
if len(prediction_paths) != len(set(prediction_paths)):
    raise ValueError("Duplicate input table paths")

votes = pd.DataFrame(index=adata.obs_names)
annotations = pd.DataFrame(index=adata.obs_names)
for method, prediction_path in zip(methods, prediction_paths):
    prediction = load(prediction_path, method)
    vote = boolean_votes(prediction.iloc[:, 0], f"{method} prediction '{prediction_path}'")
    votes[method] = vote
    annotations[f"{method}_predicted_doublet"] = vote
    annotations[f"{method}_class"] = vote.map({False: "singlet", True: "doublet"})

collisions = annotations.columns.intersection(adata.obs.columns).tolist()
if collisions:
    raise ValueError(f"Canonical annotation names collide with input observation columns: {collisions}")
annotations = annotations.reindex(sorted(annotations.columns), axis=1)

if remove_doublets:
    adata = adata[~(votes.sum(axis=1) >= threshold), :].copy()

adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")
annotations.to_pickle(f"{prefix}.obs.pkl")

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
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__,
    }
}
with open("versions.yml", "w") as handle:
    yaml.dump(versions, handle)
