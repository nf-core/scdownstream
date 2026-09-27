#!/usr/bin/env python3

import platform

import scanpy as sc
import yaml

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

adata_norm = adata.copy()
adata_norm.X = adata.layers["counts"]
sc.pp.normalize_total(adata_norm, target_sum=None)
sc.pp.log1p(adata_norm)

adata.layers["log1p"] = adata_norm.X.copy()

adata.write_h5ad(f"{prefix}.h5ad")

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
