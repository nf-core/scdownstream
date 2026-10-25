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

labels = adata.obs[column].unique()
names = {}
errors = []
for value in labels:
    if not isinstance(value, str):
        errors.append(f"{value!r} -> <unrepresentable>: labels must be non-missing strings")
        continue
    name = value.replace(" ", "_")
    if (
        name in ("", ".", "..")
        or "/" in name
        or chr(92) in name
        or any(unicodedata.category(char) == "Cc" for char in name)
    ):
        errors.append(f"{value!r} -> {name!r}: unsafe subset basename")
    else:
        destination = Path(f"{name}.h5ad")
        if destination.exists() or destination.is_symlink():
            errors.append(f"{value!r} -> {destination.name!r}: subset output path already exists")
    names.setdefault(name, []).append(value)

for name, originals in names.items():
    if len(originals) > 1:
        errors.append(f"{originals!r} -> {name!r}: subset basename collision")

if errors:
    raise ValueError(f"Invalid subset labels in column {column!r}: " + "; ".join(errors))

for name, (value,) in names.items():
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
