#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os

os.environ["KMP_AFFINITY"] = "disabled"

import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

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

kwargs = {
    "random_state": 0,
    "key_added": key_added,
}

if input_layer and input_layer in adata.layers:
    kwargs["layer"] = input_layer
else:
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)

sc.pp.pca(adata, **kwargs)

adata.write_h5ad(f"{prefix}.h5ad")
df = pd.DataFrame(adata.obsm[key_added], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions
versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
