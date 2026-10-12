#!/usr/bin/env python3

import platform

import pandas as pd
import scanpy as sc
import scimilarity
import yaml
from scimilarity import CellQuery
from scimilarity.utils import align_dataset, lognorm_counts

adata = sc.read_h5ad("${h5ad}")
adata_raw = adata.copy()

use_gpu = "${task.accelerator ? 'true' : 'false'}" == "true"
cq = CellQuery("${model}", use_gpu=use_gpu)

adata.layers["counts"] = adata.X
adata = align_dataset(adata, cq.gene_order)
adata = lognorm_counts(adata)

embeddings = cq.get_embeddings(adata.X)

# Store the embeddings
adata_raw.obsm["X_emb"] = embeddings

# Write the output
adata_raw.write_h5ad("${prefix}.h5ad")
df = pd.DataFrame(embeddings, index=adata_raw.obs_names)
df.to_pickle("X_${prefix}.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scimilarity": scimilarity.__version__,
        "scanpy": sc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
