#!/usr/bin/env python3

import platform

import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


adata = sc.read_h5ad("${h5ad}", backed="r")
prefix = "${prefix}"

kwargs = {"use_rep": "${rep}"}
n_pcs = "${n_pcs ?: ''}"
if n_pcs:
    kwargs["n_pcs"] = int(n_pcs)

sc.pp.neighbors(adata, **kwargs)

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
