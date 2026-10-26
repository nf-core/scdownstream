#!/usr/bin/env python3

import importlib.metadata
import platform

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
batch_col = "${batch_col}"
input_layer = "${input_layer}"
log_normalize = "${log_normalize}" == "true"

if not adata.obs_names.is_unique or adata.obs_names.hasnans:
    raise ValueError("Scanorama requires unique, non-missing observation identifiers.")
if batch_col not in adata.obs:
    raise ValueError(f"Scanorama batch column {batch_col!r} is not present in obs.")
if adata.obs[batch_col].isna().any():
    raise ValueError(f"Scanorama batch column {batch_col!r} contains missing values.")

if input_layer != "X" and input_layer not in adata.layers:
    raise ValueError(
        f"input_layer {input_layer!r} is not present in adata.layers (available: {list(adata.layers.keys())})"
    )

working_layer = "_scanorama_work"
if working_layer in adata.layers:
    raise ValueError(f"Scanorama working layer {working_layer!r} already exists")
source = adata.X if input_layer == "X" else adata.layers[input_layer]
adata.layers[working_layer] = source.copy()

if log_normalize:
    target_sum = float(np.median(np.asarray(adata.layers[working_layer].sum(axis=1)).ravel()))
    sc.pp.normalize_total(adata, layer=working_layer, target_sum=target_sum)
    adata.layers[working_layer] = sc.pp.log1p(adata.layers[working_layer])

adata.layers[working_layer] = sc.pp.scale(adata.layers[working_layer], max_value=10)
pca = sc.pp.pca(
    adata.layers[working_layer],
    n_comps=30,
    zero_center=False,
    mask_var=adata.var["highly_variable"].to_numpy() if "highly_variable" in adata.var else None,
)
del adata.layers[working_layer]

# Preserve first-seen batch order, including the existing contiguous-batch path.
batch_codes, _ = pd.factorize(adata.obs[batch_col], sort=False)
permutation = np.argsort(batch_codes, kind="stable")
adata_proc = sc.AnnData(pca[permutation], obs=adata.obs[[batch_col]].iloc[permutation].copy())
adata_proc.obsm["X_pca"] = adata_proc.X
sce.pp.scanorama_integrate(adata_proc, batch_col)

adata_proc = adata_proc[adata.obs_names]
embedding = adata_proc.obsm["X_scanorama"]
adata.obsm["X_scanorama"] = embedding
adata.obsm["X_emb"] = embedding
adata.write_h5ad(f"{prefix}.h5ad", compression="gzip")

pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names).to_pickle(f"X_{prefix}.pkl")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "scanorama": importlib.metadata.version("scanorama"),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
