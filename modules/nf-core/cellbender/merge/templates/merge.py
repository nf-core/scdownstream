#!/usr/bin/env python3

import platform

import anndata as ad
import cellbender
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


adata = ad.read_h5ad("${filtered}")

adata_cellbender = load_anndata_from_input_and_output("${unfiltered}", "${cellbender_h5}", analyzed_barcodes_only=False)

# Subset to the barcodes and genes present in the filtered matrix.
# Gene symbols (var index) may not be unique, so prefer Ensembl IDs when present.
# Column names differ: 10x/readh5 uses 'gene_ids'; unify/cellbender uses 'gene_id'.
filtered_gene_id_col = next((col for col in ("gene_ids", "gene_id") if col in adata.var.columns), None)
cellbender_gene_id_col = next((col for col in ("gene_id", "gene_ids") if col in adata_cellbender.var.columns), None)
filtered_ids = adata.var[filtered_gene_id_col] if filtered_gene_id_col else adata.var.index
cellbender_ids = (
    adata_cellbender.var[cellbender_gene_id_col]
    if filtered_gene_id_col and cellbender_gene_id_col
    else adata_cellbender.var.index
)
cb_id_to_pos = {gid: i for i, gid in enumerate(cellbender_ids)}
var_positions = [cb_id_to_pos[gid] for gid in filtered_ids]
adata_cellbender = adata_cellbender[adata.obs_names, var_positions]

if "${output_layer}" == "X":
    adata.X = adata_cellbender.layers["cellbender"]
else:
    adata.layers["${output_layer}"] = adata_cellbender.layers["cellbender"]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "cellbender": cellbender.__version__}}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
