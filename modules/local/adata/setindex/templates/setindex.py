#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import importlib.metadata
import anndata as ad
import yaml

adata = ad.read_h5ad("$h5ad")

axis = "$axis"
column = "$column"

df = adata.obs if axis == "obs" else adata.var

assert column in df.columns, f"Column '{column}' not found in {axis}"
df.index = df[column]
df.index.name = None

if axis == "obs":
    adata.obs = df
else:
    adata.var = df

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
