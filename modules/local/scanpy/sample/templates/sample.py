#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import yaml
import scanpy as sc
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
n = "${n ?: ''}"
fraction = "${fraction ?: ''}"

kwargs = {
    "rng": 0
}

if n:
    kwargs["n"] = int(n)
elif fraction:
    kwargs["fraction"] = float(fraction)
else:
    raise ValueError("Either n or fraction must be set")

if "n" in kwargs and kwargs["n"] >= adata.n_obs:
    print(f"Warning: n is greater than the number of cells in the dataset ({adata.n_obs}). Using the entire dataset.")
else:
    sc.pp.sample(adata, **kwargs)

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
