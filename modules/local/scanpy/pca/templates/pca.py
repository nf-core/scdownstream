#!/usr/bin/env python3

import platform

import pandas as pd
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
key_added = "${key_added}"
input_layer = "${input_layer}"
log_normalize = "${log_normalize}" == "true"

if input_layer != "X" and input_layer not in adata.layers:
    raise ValueError(
        f"input_layer {input_layer!r} is not present in adata.layers (available: {list(adata.layers.keys())})"
    )

analysis_layer = None if input_layer == "X" else input_layer
if log_normalize:
    output_layer = "pca_log1p"
    if output_layer in adata.layers:
        raise ValueError(f"Normalisation output layer {output_layer!r} already exists")
    adata.layers[output_layer] = sc.pp.normalize_total(adata, layer=analysis_layer, inplace=False)["X"]
    sc.pp.log1p(adata, layer=output_layer)
    analysis_layer = output_layer

sc.pp.pca(
    adata,
    layer=analysis_layer,
    random_state=0,
    key_added=key_added,
)

adata.write_h5ad(f"{prefix}.h5ad")
df = pd.DataFrame(adata.obsm[key_added], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions
versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
