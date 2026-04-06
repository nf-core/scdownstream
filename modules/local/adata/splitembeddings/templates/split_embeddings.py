#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform

import anndata as ad
import yaml

adata = ad.read_h5ad("${h5ad}")
embeddings = "${embeddings.join(' ')}".split()

assert len(embeddings) > 0, "No embeddings to split."
for embedding in embeddings:
    assert f"X_{embedding}" in adata.obsm, f"Embedding {embedding} not found in adata."

obsm = adata.obsm
adata.obsm = {}

for embedding in embeddings:
    adata.obsm["X_emb"] = obsm[f"X_{embedding}"]
    adata.write_h5ad(f"{embedding}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
