#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import bbknn
import pandas as pd
import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
input_layer = "${input_layer}"

pca_kwargs = {}

if input_layer and input_layer in adata.layers:
    pca_kwargs["layer"] = input_layer
else:
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)

sc.tl.pca(adata, **pca_kwargs)

kwargs = {
    "batch_key": "${batch_col}",
    "copy": True,
}

if adata.n_obs >= 1e5:
    kwargs["neighbors_within_batch"] = 25

adata = bbknn.bbknn(adata, **kwargs)

adata.write_h5ad(f"{prefix}.h5ad")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "bbknn": bbknn.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
