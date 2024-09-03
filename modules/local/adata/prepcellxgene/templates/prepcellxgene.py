#!/usr/bin/env python3

import platform
import anndata as ad
from scipy.sparse import csc_matrix
import numpy as np
import scipy as sp

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

adata = ad.read_h5ad("${h5ad}")

integration_methods = ["harmony", "scvi", "scanvi", "seurat", "bbknn", "combat"]

for integration in integration_methods:
    embedding_key = f"X_{integration}"
    if embedding_key in adata.obsm.keys():
        adata.obsm[integration] = adata.obsm.pop(embedding_key)

for layer in adata.layers.keys():
    adata.layers[layer] = csc_matrix(adata.layers[layer]).astype(np.float32)
adata.X = csc_matrix(adata.X).astype(np.float32)

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scipy": sp.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
