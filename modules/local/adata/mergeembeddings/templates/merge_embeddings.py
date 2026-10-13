#!/usr/bin/env python3

import platform

import anndata as ad
import numpy as np
import pandas as pd
import yaml


def read_anndata(path, backed=None):
    if path.rstrip("/").endswith(".zarr"):
        return ad.read_zarr(path)
    return ad.read_h5ad(path, backed=backed)


def write_anndata(adata, stem, fmt):
    if fmt == "h5ad":
        adata.write_h5ad(f"{stem}.h5ad")
    elif fmt == "zarr":
        if adata.isbacked:
            adata = adata.to_memory()
        # anndata 0.12 defaults to Zarr v2 and anndata 0.13 fails to import with ANNDATA_ZARR_WRITE_FORMAT set
        ad.settings.zarr_write_format = 3
        adata.write_zarr(f"{stem}.zarr")


adata_integrated = read_anndata("${integrated}", backed="r")
adata_base = read_anndata("${base}", backed="r")
adata_combined = read_anndata("${combined}")
integration = "${integration_key}"

emb = pd.concat(
    [
        pd.DataFrame(adata_base.obsm[f"X_{integration}"], index=adata_base.obs_names),
        pd.DataFrame(adata_integrated.obsm["X_emb"], index=adata_integrated.obs_names),
    ],
    axis=0,
)

adata_combined.obsm["X_emb"] = emb.loc[adata_combined.obs_names].to_numpy()

# Drop all rows with NaN in the embeddings
# These only occur when during the creation of the base anndata there were cells with no expression in the highly variable genes
# Needs to be done because downstream processes (e.g. clustering) require the embeddings to be complete
adata_combined = adata_combined[~np.isnan(adata_combined.obsm["X_emb"]).any(axis=1)]

if integration == "scanvi":
    labels = pd.concat(
        [
            pd.DataFrame(adata_base.obs["label:scANVI"], index=adata_base.obs_names),
            pd.DataFrame(adata_integrated.obs["label:scANVI"], index=adata_integrated.obs_names),
        ],
        axis=0,
    )

    adata_combined.obs["label:scANVI"] = labels.loc[adata_combined.obs_names]

    adata_combined.obs[["label:scANVI"]].to_pickle("${prefix}.pkl")

df = pd.DataFrame(adata_combined.obsm["X_emb"], index=adata_combined.obs_names)
df.to_pickle("X_${prefix}.pkl")

write_anndata(adata_combined, "${prefix}", "${anndata_format}")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "anndata": ad.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
