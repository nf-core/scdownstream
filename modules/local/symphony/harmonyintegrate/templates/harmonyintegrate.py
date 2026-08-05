#!/usr/bin/env python3

import os

os.environ["KMP_AFFINITY"] = "disabled"
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import importlib.metadata
import platform

import numpy as np
import pandas as pd
import scanpy as sc
import symphonypy as sp
import yaml
from anndata import AnnData
from scipy.sparse import csr_matrix
from threadpoolctl import threadpool_limits


def build_reference(adata, target_sum):
    harmony = adata.uns["harmony"]
    return AnnData(
        X=csr_matrix((0, adata.n_vars), dtype=np.float32),
        var=adata.var[["mean", "std", "highly_variable"]].copy(),
        varm={"PCs": adata.varm["PCs"].copy()},
        uns={
            "harmony": {
                "Nr": harmony["Nr"],
                "C": harmony["C"],
                "K": harmony["K"],
                "sigma": harmony.get("sigma"),
                "ref_basis_loadings": harmony["ref_basis_loadings"],
            },
            "normalize": {"target_sum": target_sum},
        },
    )


threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
adata_proc = adata.copy()
prefix = "${prefix}"
batch_col = "${batch_col}"
counts_layer = "${counts_layer}"

if counts_layer != "X":
    adata_proc.X = adata_proc.layers[counts_layer]

target_sum = float(np.median(np.asarray(adata_proc.X.sum(axis=1)).ravel()))
sc.pp.normalize_total(adata_proc, target_sum=target_sum)
sc.pp.log1p(adata_proc)
sc.pp.scale(adata_proc, max_value=10)
sc.pp.pca(adata_proc, n_comps=30, zero_center=False)
if "highly_variable" not in adata_proc.var.columns:
    adata_proc.var["highly_variable"] = True

sp.pp.harmony_integrate(
    adata_proc,
    key=batch_col,
    flavor="python",
    ref_basis_source="X_pca",
    ref_basis_adjusted="X_symphony",
)

adata_proc.uns["symphony"] = adata_proc.uns["harmony"]
adata_proc.uns["normalize"] = {"target_sum": target_sum}

build_reference(adata_proc, target_sum).write_h5ad(f"{prefix}_reference.h5ad")

adata.obsm["X_symphony"] = adata_proc.obsm["X_symphony"]
adata.obsm["X_emb"] = adata_proc.obsm["X_symphony"]
adata.write_h5ad(f"{prefix}.h5ad")

pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names).to_parquet(f"X_{prefix}.parquet", index=True)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": importlib.metadata.version("scanpy"),
        "symphonypy": importlib.metadata.version("symphonypy"),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
