#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

import platform
import argparse
import shlex

import numpy as np

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


adata = sc.read_h5ad("${h5ad}", backed='r')
prefix = "${prefix}"
args = "${args}"
parser = argparse.ArgumentParser()
parser.add_argument("--decimals", type=int, default=None)
params = parser.parse_args(shlex.split(args))

kwargs = {
    "use_rep": "${rep}"
}

sc.pp.neighbors(adata, **kwargs)

if params.decimals is not None:
    for key in adata.obsp:
        if hasattr(adata.obsp[key], 'data'):
            adata.obsp[key].data = np.round(adata.obsp[key].data.astype(np.float64), params.decimals)

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
