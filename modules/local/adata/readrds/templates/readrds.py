#!/usr/bin/env python3

import os

os.environ["TMPDIR"] = "."

import anndata as ad
import anndata2ri
import rpy2
import pandas as pd
import rpy2.robjects as ro
seurat = ro.packages.importr('Seurat')
# Import SingleCellExperiment to check for class
sce_pkg = ro.packages.importr('SingleCellExperiment')

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

# Read the RDS file first
rds_obj = ro.r(f'readRDS("${rds}")')

# Check if it's already a SingleCellExperiment
is_sce = 'SingleCellExperiment' in rds_obj.extends()

# Convert only if not already a SingleCellExperiment
if not is_sce:
    sce = seurat.as_SingleCellExperiment(rds_obj)
else:
    sce = rds_obj

adata = anndata2ri.rpy2py(sce)

# Convert indices to string
adata.obs.index = adata.obs.index.astype(str)
adata.var.index = adata.var.index.astype(str)

adata.write_h5ad("${prefix}.h5ad")

versions = {
    "${task.process}": {
        "anndata": ad.__version__,
        "anndata2ri": anndata2ri.__version__,
        "rpy2": rpy2.__version__,
        "pandas": pd.__version__,
        "seurat": seurat.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
