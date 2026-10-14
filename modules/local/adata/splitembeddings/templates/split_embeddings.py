#!/usr/bin/env python3

import platform
from pathlib import Path

import anndata as ad
import yaml


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        # An absolute path keeps zarr from parsing colons in Nextflow prefixes as URL schemes
        return ad.read_zarr(Path(path).absolute())
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(Path(f"{stem}.zarr").absolute())


adata = read_anndata("${anndata}")
embeddings = "${embeddings.join(' ')}".split()

assert len(embeddings) > 0, "No embeddings to split."
for embedding in embeddings:
    assert f"X_{embedding}" in adata.obsm, f"Embedding {embedding} not found in adata."

obsm = adata.obsm
adata.obsm = {}

for embedding in embeddings:
    adata.obsm["X_emb"] = obsm[f"X_{embedding}"]
    write_anndata(adata, embedding, "${anndata_format}")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "anndata": ad.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
