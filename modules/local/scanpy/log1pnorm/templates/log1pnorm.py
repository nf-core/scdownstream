#!/usr/bin/env python3

import platform

import scanpy as sc
import yaml

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

adata.layers["log1p"] = adata.X.copy()
sc.pp.normalize_total(adata, layer="log1p", target_sum=None)
sc.pp.log1p(adata, layer="log1p")

adata.write_h5ad(f"{prefix}.h5ad")

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
