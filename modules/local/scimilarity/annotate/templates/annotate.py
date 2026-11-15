#!/usr/bin/env python3

import platform
from pathlib import Path

import pyarrow
import scanpy as sc
import scimilarity
import yaml
from scimilarity import CellAnnotation

adata = sc.read_h5ad("${h5ad}")

use_gpu = "${task.accelerator ? 'true' : 'false'}" == "true"
ca = CellAnnotation("${model}", use_gpu=use_gpu)

predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(adata.obsm["X_emb"])

adata.obs["annotation:scimilarity"] = predictions.values

# Write the output
adata.write_h5ad("${prefix}.h5ad", compression="gzip")
obs_fragment = adata.obs[["annotation:scimilarity"]].copy()
if "__index__" in obs_fragment.columns or not obs_fragment.columns.is_unique:
    raise ValueError("Observation fragment columns must be unique and must not contain reserved '__index__'")
if (
    obs_fragment.index.hasnans
    or not obs_fragment.index.is_unique
    or not all(isinstance(value, str) for value in obs_fragment.index)
):
    raise ValueError("Observation fragment index must contain unique, non-null strings")
obs_fragment.reset_index(names="__index__").to_parquet(
    Path("${prefix}.parquet").resolve(), engine="pyarrow", index=False, compression="snappy"
)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scimilarity": scimilarity.__version__,
        "scanpy": sc.__version__,
        "pyarrow": pyarrow.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
