#!/usr/bin/env python3

import os
import shutil
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"

import scanpy as sc
from scimilarity import CellAnnotation
import scimilarity
import yaml

adata = sc.read_h5ad("${h5ad}")
adata_raw = adata.copy()

use_gpu = "${task.ext.use_gpu}" == "true"
ca = CellAnnotation("${model}", use_gpu=use_gpu)

predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(
    adata.obsm["X_emb"]
)

adata_raw.obs["annotation:scimilarity"] = predictions.values

# Write the output
adata_raw.write_h5ad("${prefix}.h5ad")
adata_raw.obs[["annotation:scimilarity"]].to_pickle("${prefix}.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scimilarity": scimilarity.__version__,
        "scanpy": sc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
