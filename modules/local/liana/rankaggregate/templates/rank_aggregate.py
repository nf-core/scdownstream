#!/usr/bin/env python3

import os
import platform
import yaml

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import pandas as pd
import scanpy as sc
import liana as li

from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
obs_key = "${obs_key}"

if adata.obs[obs_key].nunique() > 1:
    if (adata.X < 0).nnz == 0:
        sc.pp.log1p(adata)
    try:
        li.mt.rank_aggregate(
            adata, obs_key, use_raw=False, verbose=True, n_jobs=int("${task.cpus}")
        )
        df: pd.DataFrame = adata.uns["liana_res"]

        df.to_pickle(f"{prefix}.pkl")
        adata.write_h5ad(f"{prefix}.h5ad")

    except ValueError as e:
        if "cannot set a frame with no defined index and a scalar" in str(e):
            print(f"Error: {e}")
        else:
            raise e
else:
    print(
        f"Skipping rank aggregation because the column {obs_key} has only one unique value."
    )

# Versions

versions = {
    "python": platform.python_version(),
    "scanpy": sc.__version__,
    "liana": li.__version__,
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
