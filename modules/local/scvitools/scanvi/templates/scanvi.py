#!/usr/bin/env python3

import os

os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

import scvi
import anndata as ad
import pandas as pd
from scvi.model import SCVI, SCANVI
import platform
import torch
import yaml

torch.use_deterministic_algorithms(True)
torch.set_float32_matmul_precision('medium')

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

scvi.settings.num_threads = int("${task.cpus}")
scvi.settings.seed = 0

adata = ad.read_h5ad("${h5ad}")
reference_model_path = "reference_model"
reference_model_type = "${meta2.id}"

# FIXME: This is a hack, the columns from celltypist are:
# f"celltypist:{model_name}", f"celltypist:{model_name}:conf" (from line 69)
# But here it expects the labels under "label" currently it seems that no
# script triggered is handling confidence cutoff and label assignment to
# "${label_col}" ("label").
# Here I hack a solution in the event celltypist is used:
# Overwrite labels with CellTypist predictions (with confidence filter)
# Find any CellTypist prediction/confidence columns
pred_cols = [
    c
    for c in adata.obs.columns
    if c.startswith("celltypist:") and not c.endswith(":conf")
]

# If multiple models are present, I will pick the first one
pred_col = pred_cols[0]
conf_col = pred_col + ":conf"

print(f"Using CellTypist column: {pred_col}")

if conf_col in adata.obs:
    CUTOFF = 0.7
    adata.obs["${label_col}"] = adata.obs.apply(
        lambda r: r[pred_col] if r[conf_col] >= CUTOFF else "unknown", axis=1
    )
else:
    adata.obs["${label_col}"] = adata.obs[pred_col]

# Make sure labels are categorical and contain "unknown"
adata.obs["${label_col}"] = adata.obs["${label_col}"].astype("category")
if "unknown" not in adata.obs["${label_col}"].cat.categories:
    adata.obs["${label_col}"] = adata.obs["${label_col}"].cat.add_categories(["unknown"])
# ---- End Hack ----

if reference_model_type == "scanvi":
    SCANVI.prepare_query_anndata(adata, reference_model_path)
    model = SCANVI.load_query_data(adata, reference_model_path)
else:
    unique_labels = set(adata.obs["${label_col}"].unique())
    unique_labels.discard("unknown")

    if not len(unique_labels) > 1:
        raise ValueError("Not enough labels to run scANVI.")

    if reference_model_type == "scvi":
        SCVI.prepare_query_anndata(adata, reference_model_path)
        model = SCVI.load(reference_model_path, adata)
        model = SCANVI.from_scvi_model(
            scvi_model=model, labels_key="${label_col}", unlabeled_category="unknown"
        )
    else:
        categorical_covariates = "${categorical_covariates}"
        continuous_covariates = "${continuous_covariates}"

        categorical_covariates = categorical_covariates.split(",") if categorical_covariates else None
        continuous_covariates = continuous_covariates.split(",") if continuous_covariates else None

        SCANVI.setup_anndata(adata, batch_key="${batch_col}", labels_key="${label_col}", unlabeled_category="${unlabeled_category}",
                                categorical_covariate_keys = categorical_covariates,
                                continuous_covariate_keys = continuous_covariates)

        model = SCANVI(adata,
                        n_hidden=int("${n_hidden}"),
                        n_layers=int("${n_layers}"),
                        n_latent=int("${n_latent}"),
                        dispersion="${dispersion}",
                        gene_likelihood="${gene_likelihood}")

if "${task.ext.use_gpu}" == "true":
    model.to_device(0)

model.train(early_stopping=True,
            max_epochs=int("${max_epochs}") if "${max_epochs?:''}" else None)

# Round to ensure hashes are stable
adata.obsm["X_emb"] = model.get_latent_representation()

adata.obs["label:scANVI"] = model.predict()

del adata.uns["_scvi_manager_uuid"]
del adata.uns["_scvi_uuid"]

adata.write_h5ad("${prefix}.h5ad")
adata.obs[["label:scANVI"]].to_pickle("${prefix}.pkl")
model.save("${prefix}_model")

df = pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names)
df.to_pickle("X_${prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "scvi": scvi.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
