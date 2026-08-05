#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["CELLTYPIST_FOLDER"] = "./tmp/celltypist"

import celltypist
import pandas as pd
import scanpy as sc
from celltypist import models as ct_models


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data: The dictionary to format.
        indent: The current indentation level.

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


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

models = "${models.join(' ')}".split()

adata_celltypist = adata.copy()  # make a copy of our adata
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

symbol_col = "${symbol_col}"
if symbol_col != "index" and symbol_col:
    if symbol_col not in adata_celltypist.var.columns:
        raise ValueError(f"Symbol column {symbol_col} not found in adata.var.columns")
    adata_celltypist.var_names = adata_celltypist.var[symbol_col]

# celltypist expects a string index, because it will make unique names by appending "-1", "-2"
# to duplicates if necessary. Cast other types (e.g. CategoricalIndex) to str:
adata_celltypist.var_names = adata_celltypist.var_names.astype(str)

df_list = []
manifest_rows = []

for model in models:
    model_file = f"{model}.pkl" if not model.endswith(".pkl") else model
    model_name = model_file[:-4]
    ct_models.download_models(model=model_file)
    model_obj = ct_models.Model.load(model_file)

    predictions = celltypist.annotate(adata_celltypist, model=model_obj)
    predictions_adata = predictions.to_adata()

    per_cell_col = f"annotation:celltypist:{model_name}:per_cell"
    conf_col = f"annotation:celltypist:{model_name}:per_cell:confidence"
    df_celltypist = predictions_adata.obs.loc[adata.obs.index, ["predicted_labels", "conf_score"]]
    df_celltypist.columns = [per_cell_col, conf_col]
    df_list.append(df_celltypist)
    manifest_rows.extend(
        [
            {"obs_column": per_cell_col, "aggregatable": "true"},
            {"obs_column": conf_col, "aggregatable": "false"},
        ]
    )

df_celltypist = pd.concat(df_list, axis=1)
df_celltypist.to_parquet("${prefix}.parquet", index=True)

pd.DataFrame(manifest_rows).to_csv(f"{prefix}_annotation_columns.csv", index=False)

adata.obs = pd.concat([adata.obs, df_celltypist], axis=1)
adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "celltypist": celltypist.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
