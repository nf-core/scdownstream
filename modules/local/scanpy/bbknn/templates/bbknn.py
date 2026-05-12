#!/usr/bin/env python3

import os
import platform
import argparse
import shlex

import numpy as np

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import pandas as pd
import bbknn
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")

sc.tl.pca(adata)

kwargs = {
    "batch_key": "${batch_col}",
    "copy": True,
}

if adata.n_obs >= 1e5:
    kwargs["neighbors_within_batch"] = 25

adata = bbknn.bbknn(adata, **kwargs)

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "bbknn": bbknn.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
