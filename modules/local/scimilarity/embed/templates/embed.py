#!/usr/bin/env python3

import platform
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scimilarity
import yaml
from scimilarity import CellQuery
from scimilarity.utils import align_dataset

adata = sc.read_h5ad("${h5ad}")

use_gpu = "${task.accelerator ? 'true' : 'false'}" == "true"
cq = CellQuery("${model}", use_gpu=use_gpu)

aligned = align_dataset(
    ad.AnnData(X=adata.X, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names)),
    cq.gene_order,
)
aligned.layers["log1p"] = aligned.X.copy()
if isinstance(aligned.layers["log1p"], np.ndarray) and np.isnan(aligned.layers["log1p"]).any():
    warnings.warn("NaNs detected in expression counts; replacing them with zero before normalisation.", UserWarning)
    aligned.layers["log1p"] = np.nan_to_num(aligned.layers["log1p"], nan=0.0)
sc.pp.normalize_total(aligned, layer="log1p", target_sum=1e4)
sc.pp.log1p(aligned, layer="log1p")

embeddings = cq.get_embeddings(aligned.layers["log1p"])

# Store the embeddings
adata.obsm["X_emb"] = embeddings

# Write the output
adata.write_h5ad("${prefix}.h5ad")
df = pd.DataFrame(embeddings, index=adata.obs_names)
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
