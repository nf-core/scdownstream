#!/usr/bin/env python3

import scanpy as sc
import scipy
from scipy.sparse import csr_matrix
import platform

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

adata = sc.read_h5ad("$h5ad")

# Unify batches
batch_col = "${meta.batch_col}"
if batch_col not in adata.obs:
    adata.obs[batch_col] = "${meta.id}"

if batch_col != "batch":
    if "batch" in adata.obs:
        raise ValueError("The batch column already exists.")
    adata.obs["batch"] = adata.obs[batch_col]
    del adata.obs[batch_col]

# Unify labels
label_col = "${meta.label_col ?: ''}"
unknown_label = "${meta.unknown_label}"

if label_col:
    if label_col != "label":
        if "label" in adata.obs:
            raise ValueError("The label column already exists.")
        adata.obs["label"] = adata.obs[label_col]
        del adata.obs[label_col]

    if unknown_label != "unknown":
        if "unknown" in adata.obs["label"]:
            raise ValueError("The label column already contains 'unknown' values.")
        adata.obs["label"].replace({unknown_label: "unknown"}, inplace=True)
else:
    if "label" in adata.obs:
        raise ValueError("The label column already exists.")
    adata.obs["label"] = "unknown"

# Add "sample" column
if "sample" in adata.obs and not adata.obs["sample"].equals("${meta.id}"):
    adata.obs["sample_original"] = adata.obs["sample"]
adata.obs["sample"] = "${meta.id}"

# Convert to CSR matrix
adata.X = csr_matrix(adata.X)
adata.layers["counts"] = adata.X

# Perform basic filtering
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "scipy": scipy.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
