#!/usr/bin/env python3

import scvi
import anndata as ad
from scvi.model import SCVI
from scvi.external import SOLO

adata = ad.read_h5ad("${h5ad}")

setup_kwargs = {"layer": "counts", "batch_key": "${meta.batch_col}"}

# Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
model_kwargs = {
    "gene_likelihood": "nb",
    "n_layers": 2,
    "n_hidden": 128,
    "n_latent": 30,
}

train_kwargs = {"train_size": 1.0}

SCVI.setup_anndata(adata, **setup_kwargs)
model = SCVI(adata, **model_kwargs)
model.train(**train_kwargs)

if "${meta.label_col}":
    from scvi.model import SCANVI

    model = SCANVI.from_scvi_model(scvi_model = model, labels_key = "${meta.label_col}", unlabeled_category = "Unknown")
    model.train(**train_kwargs)

adata.obs["singlet"] = False

for batch in adata.obs["${meta.batch_col}"].unique():
    solo = SOLO.from_scvi_model(model, restrict_to_batch=batch)
    solo.train()
    result = solo.predict(False)

    singlets = result[result == "singlet"].index.tolist()
    adata.obs.loc[singlets, "singlet"] = True

adata = adata[adata.obs["singlet"]].copy()
adata.obs.drop("singlet", axis=1, inplace=True)

adata.write_h5ad("${prefix}.h5ad")
