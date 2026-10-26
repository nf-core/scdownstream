#!/usr/bin/env python3

import importlib.metadata
import platform
import unicodedata
from pathlib import Path

import anndata as ad
import yaml

adata = ad.read_h5ad("${h5ad}")
column = "${column}"

assert column in adata.obs.columns, f"Column {column} not found in adata."

names = {}
for value in adata.obs[column].unique():
    if not isinstance(value, str):
        raise ValueError(f"{value!r}: subset labels must be non-missing strings")
    name = value.replace(" ", "_")
    if (
        name in ("", ".", "..")
        or "/" in name
        or chr(92) in name  # Backslash
        or any(unicodedata.category(char) == "Cc" for char in name)
    ):
        raise ValueError(f"{value!r} -> {name!r}: unsafe subset basename")
    if name in names:
        raise ValueError(f"{[names[name], value]!r} -> {name!r}: subset basename collision")
    destination = Path(f"{name}.h5ad")
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"{value!r} -> {destination.name!r}: subset output path already exists")
    names[name] = value

for name, value in names.items():
    adata_subset = adata[adata.obs[column] == value].copy()
    adata_subset.write_h5ad(f"{name}.h5ad", compression="gzip")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
