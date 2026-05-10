#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"

import pandas as pd
import scanpy as sc
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellQuery
import scimilarity
import yaml

adata = sc.read_h5ad("${h5ad}")
adata_raw = adata.copy()

use_gpu = "${task.ext.use_gpu}" == "true"
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
df.to_pickle("X_${meta.id}.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scimilarity": scimilarity.__version__,
        "scanpy": sc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
