#!/usr/bin/env python3

import platform

import anndata as ad
import cellbender
import numpy as np
import pandas as pd
from cellbender.remove_background.downstream import load_anndata_from_input_and_output


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


def var_alignment_keys(adata):
    if "probe_ids" in adata.var.columns:
        keys = adata.var["probe_ids"].astype(str)
        if keys.is_unique:
            return keys
    if "gene_ids" in adata.var.columns:
        keys = adata.var["gene_ids"].astype(str)
        if keys.is_unique:
            return keys
    keys = adata.var.index.astype(str)
    if keys.is_unique:
        return keys
    raise ValueError(
        "Cannot align CellBender output to the filtered matrix: no unique feature "
        "identifiers found among var `probe_ids`, `gene_ids`, or the var index."
    )


def subset_cellbender_to_filtered(adata, adata_cellbender):
    adata_cellbender = adata_cellbender[adata.obs_names]
    if adata.n_vars == adata_cellbender.n_vars:
        return adata_cellbender

    target_keys = var_alignment_keys(adata)
    source_keys = var_alignment_keys(adata_cellbender)

    source_positions = pd.Series(np.arange(adata_cellbender.n_vars), index=source_keys)
    if not source_positions.index.is_unique:
        source_positions = source_positions[~source_positions.index.duplicated(keep="first")]

    positions = source_positions.loc[target_keys.values]
    if positions.isna().any():
        missing = target_keys[positions.isna()].unique()[:5]
        raise ValueError(
            "Features from the filtered matrix were not found in the CellBender output. "
            f"Examples: {list(missing)}"
        )

    return adata_cellbender[:, positions.to_numpy(dtype=int)]


adata = ad.read_h5ad("${filtered}")

adata_cellbender = load_anndata_from_input_and_output("${unfiltered}", "${cellbender_h5}", analyzed_barcodes_only=False)

adata_cellbender = subset_cellbender_to_filtered(adata, adata_cellbender)

if "${output_layer}" == "X":
    adata.X = adata_cellbender.layers["cellbender"]
else:
    adata.layers["${output_layer}"] = adata_cellbender.layers["cellbender"]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "cellbender": cellbender.__version__}}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
