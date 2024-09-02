#!/usr/bin/env python3

import scvi
import anndata as ad
import pandas as pd
from scvi.model.base import ArchesMixin as sca
import platform

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

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

adata_query     = ad.read_h5ad("${h5ad_transfer}")
adata_inner     = ad.read_h5ad("${h5ad_inner}")
model_reference = "model"

sca.prepare_query_anndata(
    adata_query,
    model_reference,
    return_reference_var_names=True
)

model = sca.load_query_data(
    adata_query,
    model_reference
)

if "${task.ext.use_gpu}" == "true":
    model.to_device(0)

model.train()
model.save("${prefix}_model")

embedding = model.get_latent_representation()
adata_inner.obsm.loc[adata_query.obs_names, "X_emb"] = embedding

adata_inner.write_h5ad("${prefix}.h5ad")
adata_inner.obsm[["X_emb"]].to_pickle("${prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scvi": scvi.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
