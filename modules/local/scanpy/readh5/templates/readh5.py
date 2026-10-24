#!/usr/bin/env python3

import platform

import pandas as pd
import scanpy as sc
import yaml

adata = sc.read_10x_h5("${h5}")
adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
