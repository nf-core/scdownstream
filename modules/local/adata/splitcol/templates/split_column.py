#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import importlib.metadata

import anndata as ad
import yaml

adata = ad.read_h5ad("${h5ad}")
column = "${column}"

assert column in adata.obs.columns, f"Column {column} not found in adata."

for value in adata.obs[column].unique():
    adata_subset = adata[adata.obs[column] == value]
    value = value.replace(" ", "_")
    adata_subset.write_h5ad(f"{value}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
