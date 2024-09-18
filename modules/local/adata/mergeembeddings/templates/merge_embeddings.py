#!/usr/bin/env python3

import platform

import pandas as pd
import anndata as ad

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

adata_integrated = ad.read_h5ad("${integrated}", backed="r")
adata_base = ad.read_h5ad("${base}", backed="r")
adata_combined = ad.read_h5ad("${combined}")
integration = "${meta.id}"

emb = pd.concat([
    pd.DataFrame(adata_base.obsm[f"X_{integration}"], index=adata_base.obs_names), 
    pd.DataFrame(adata_integrated.obsm["X_emb"], index=adata_integrated.obs_names)
], axis=0)

adata_combined.obsm["X_emb"] = emb.loc[adata_combined.obs_names].to_numpy()

if integration == "scanvi":
    labels = pd.concat([
        pd.DataFrame(adata_base.obs["label:scANVI"], index=adata_base.obs_names),
        pd.DataFrame(adata_integrated.obs["label:scANVI"], index=adata_integrated.obs_names)
    ], axis=0)

    adata_combined.obs["label:scANVI"] = labels.loc[adata_combined.obs_names]

    adata_combined.obs[["label:scANVI"]].to_pickle("${prefix}.pkl")

df = pd.DataFrame(adata_combined.obsm["X_emb"], index=adata_combined.obs_names)
df.to_pickle("X_${prefix}.pkl")

adata_combined.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
