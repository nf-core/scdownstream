#!/usr/bin/env python3

import os

os.environ["TMPDIR"] = "."

import importlib.metadata

import anndata as ad
import anndata2ri
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter

seurat = ro.packages.importr("Seurat")
# Import SingleCellExperiment to check for class
sce_pkg = ro.packages.importr("SingleCellExperiment")


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
rds_obj = ro.r('readRDS("${rds}")')

# Check if it's already a SingleCellExperiment
is_sce = "SingleCellExperiment" in rds_obj.extends()

# Convert only if not already a SingleCellExperiment
if not is_sce:
    sce = seurat.as_SingleCellExperiment(rds_obj)
else:
    sce = rds_obj

with localconverter(anndata2ri.converter):
    adata = ro.conversion.rpy2py(sce)

# Convert indices to string
adata.obs.index = adata.obs.index.astype(str)
adata.var.index = adata.var.index.astype(str)


def _dataframe_for_h5ad(df: pd.DataFrame) -> pd.DataFrame:
    """Rewrite string indexes/columns so nft-anndata and R can read the written H5AD."""
    df = df.copy()
    df.index = pd.CategoricalIndex(df.index.astype(str).to_list(), name=df.index.name)
    for col in df.columns:
        if isinstance(df[col].dtype, pd.StringDtype) or df[col].dtype == object:
            df[col] = pd.Series(df[col].astype(str).to_list(), index=df.index, name=col, dtype=object)
    return df


def _prepare_adata_for_h5ad(adata_obj):
    """Avoid nullable-string-array encodings that downstream tools cannot read."""
    adata_obj.obs = _dataframe_for_h5ad(adata_obj.obs)
    adata_obj.var = _dataframe_for_h5ad(adata_obj.var)
    for uns_key, uns_value in list(adata_obj.uns.items()):
        if isinstance(uns_value, pd.DataFrame):
            adata_obj.uns[uns_key] = _dataframe_for_h5ad(uns_value)
        elif isinstance(uns_value, dict):
            for key, value in list(uns_value.items()):
                if isinstance(value, pd.DataFrame):
                    uns_value[key] = _dataframe_for_h5ad(value)


_prepare_adata_for_h5ad(adata)
adata.write_h5ad("${prefix}.h5ad")

versions = {
    "${task.process}": {
        "anndata": ad.__version__,
        "anndata2ri": importlib.metadata.version("anndata2ri"),
        "rpy2": importlib.metadata.version("rpy2"),
        "pandas": pd.__version__,
        "seurat": seurat.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
