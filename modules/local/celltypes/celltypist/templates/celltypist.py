#!/usr/bin/env python3

import scanpy as sc
import celltypist
from celltypist import models
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


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
model = "${model}"

adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(
    adata_celltypist, counts_per_cell_after=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

model = f"{model}.pkl" if not model.endswith(".pkl") else model
models.download_models(model=model)
model = models.Model.load(model)

predictions = celltypist.annotate(
    adata_celltypist, model=model
)
predictions_adata = predictions.to_adata()

df_celltypist = predictions_adata.obs.loc[
    adata.obs.index, ["predicted_labels", "conf_score"]
]

df_celltypist.columns = ["cell_type:celltypist", "celltypist_confidence"]
df_celltypist.to_pickle("${prefix}.pkl")

predictions_adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "celltypist": celltypist.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
