#!/usr/bin/env python3

import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import yaml

adata_integrated = ad.read_h5ad("${integrated}", backed="r")
adata_base = ad.read_h5ad("${base}", backed="r")
adata_combined = ad.read_h5ad("${combined}")
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

    adata_combined.obs[["label:scANVI"]].reset_index(names="__index__").to_parquet(
        Path("${prefix}.parquet").resolve(), engine="pyarrow", index=False, compression="snappy"
    )

df = pd.DataFrame(adata_combined.obsm["X_emb"], index=adata_combined.obs_names)
df.columns = df.columns.astype(str)
df.reset_index(names="__index__").to_parquet(
    Path("X_${prefix}.parquet").resolve(), engine="pyarrow", index=False, compression="snappy"
)

adata_combined.write_h5ad("${prefix}.h5ad", compression="gzip")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "anndata": ad.__version__, "pyarrow": pa.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
