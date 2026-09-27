#!/usr/bin/env python3

import os
import platform
from importlib.metadata import version

# panhumanpy caches the model weights in Path.home(), which is resolved at import time
os.environ["HOME"] = os.path.abspath("./tmp")

import anndata as ad
import pandas as pd
import panhumanpy as ph
import tensorflow as tf

tf.config.threading.set_intra_op_parallelism_threads(int("${task.cpus}"))
tf.config.threading.set_inter_op_parallelism_threads(int("${task.cpus}"))


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


adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"
counts_layer = "${counts_layer}"

if counts_layer != "X" and counts_layer not in adata.layers:
    raise ValueError(f"Counts layer {counts_layer} not found in adata.layers")
counts = adata.X if counts_layer == "X" else adata.layers[counts_layer]

var = pd.DataFrame(index=adata.var_names.astype(str))
feature_names_col = None
if symbol_col != "index" and symbol_col:
    if symbol_col not in adata.var.columns:
        raise ValueError(f"Symbol column {symbol_col} not found in adata.var.columns")
    var[symbol_col] = adata.var[symbol_col].astype(str).values
    feature_names_col = symbol_col

# A minimal copy keeps panhumanpy from writing its columns into the input obs.
# panhumanpy normalises the counts itself when they are integers.
adata_query = ad.AnnData(X=counts, obs=pd.DataFrame(index=adata.obs_names), var=var)

azimuth = ph.AzimuthNN(adata_query, feature_names_col=feature_names_col, model_version="v1")
embedding = azimuth.azimuth_embed()
cells_meta = azimuth.cells_meta.loc[adata.obs_names]

columns = {
    "azimuth_broad": ("annotation:azimuth:broad:per_cell", "true"),
    "azimuth_medium": ("annotation:azimuth:medium:per_cell", "true"),
    "azimuth_fine": ("annotation:azimuth:fine:per_cell", "true"),
    "final_level_confidence": ("annotation:azimuth:fine:per_cell:confidence", "false"),
    "full_hierarchical_labels": ("annotation:azimuth:full_hierarchy:per_cell", "false"),
}

df_azimuth = cells_meta[list(columns)].rename(columns={src: dst for src, (dst, _) in columns.items()})
df_azimuth.to_pickle(f"{prefix}.pkl")

pd.DataFrame([{"obs_column": dst, "aggregatable": aggregatable} for dst, aggregatable in columns.values()]).to_csv(
    f"{prefix}_annotation_columns.csv", index=False
)

df_embedding = pd.DataFrame(embedding, index=adata.obs_names)
df_embedding.to_pickle("X_azimuth.pkl")

adata.obs = pd.concat([adata.obs, df_azimuth], axis=1)
adata.obsm["X_azimuth"] = embedding
adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "anndata": version("anndata"),
        "tensorflow": tf.__version__,
        "panhumanpy": ph.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
