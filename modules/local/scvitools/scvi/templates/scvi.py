#!/usr/bin/env python3

# Disable OpenMP CPU topology detection for MacOS compatibility
import os
os.environ["KMP_AFFINITY"] = "disabled"

os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

import scvi
import anndata as ad
import pandas as pd
from scvi.model import SCVI
import torch
import yaml

torch.set_float32_matmul_precision("medium")
torch.use_deterministic_algorithms(True)

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

scvi.settings.num_threads = int("${task.cpus}")
scvi.settings.seed = 0

adata = ad.read_h5ad("${h5ad}")
reference_model_path = "reference_model"
reference_model_type = "${meta2.id ?: ''}"

plan_kwargs = {}

if reference_model_type:
    if reference_model_type == "scanvi":
        raise ValueError("scVI does not support scANVI models.")
    elif reference_model_type == "scvi":
        SCVI.prepare_query_anndata(adata, reference_model_path)
        model = SCVI.load_query_data(adata, reference_model_path)
        plan_kwargs['weight_decay'] = 0.0
    else:
        raise ValueError(f"Invalid reference model type: {reference_model_type}")
else:
    categorical_covariates = "${categorical_covariates}"
    continuous_covariates = "${continuous_covariates}"

    categorical_covariates = (
        categorical_covariates.split(",") if categorical_covariates else None
    )
    continuous_covariates = (
        continuous_covariates.split(",") if continuous_covariates else None
    )

    SCVI.setup_anndata(
        adata,
        batch_key="${batch_col}",
        categorical_covariate_keys=categorical_covariates,
        continuous_covariate_keys=continuous_covariates,
    )

    model = SCVI(
        adata,
        n_hidden=int("${n_hidden}"),
        n_layers=int("${n_layers}"),
        n_latent=int("${n_latent}"),
        dispersion="${dispersion}",
        gene_likelihood="${gene_likelihood}",
    )

if "${task.ext.use_gpu}" == "true":
    model.to_device(0)

model.train(
    early_stopping=True,
    max_epochs=int("${max_epochs}") if "${max_epochs?:''}" else None,
    plan_kwargs=plan_kwargs,
)

# Round to ensure hashes are stable
adata.obsm["X_emb"] = model.get_latent_representation()

del adata.uns["_scvi_manager_uuid"]
del adata.uns["_scvi_uuid"]

adata.write_h5ad("${prefix}.h5ad")
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
