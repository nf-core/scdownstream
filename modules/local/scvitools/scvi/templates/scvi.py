#!/usr/bin/env python3

import os

os.environ.setdefault("TORCHINDUCTOR_CACHE_DIR", os.path.join(os.getcwd(), "torch_cache"))

os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

import anndata as ad
import pandas as pd
import scvi
import torch
import yaml
from scvi.model import SCVI
from threadpoolctl import threadpool_limits


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


torch.set_float32_matmul_precision("medium")
torch.use_deterministic_algorithms(True)

threadpool_limits(int("${task.cpus}"))

scvi.settings.num_threads = int("${task.cpus}")
scvi.settings.seed = 0

adata = read_anndata("${anndata}")
adata_work = adata.copy()
reference_model_dir = "reference_model"
reference_model_file = f"{reference_model_dir}/model.pt"
reference_model_type = "${meta2.id ?: ''}"

plan_kwargs = {}

if reference_model_type:
    if reference_model_type == "scanvi":
        raise ValueError("scVI does not support scANVI models.")
    elif reference_model_type == "scvi":
        state = torch.load(reference_model_file, map_location="cpu", weights_only=False)
        setup = state["attr_dict"]["registry_"]["setup_args"]
        batch_key = setup["batch_key"]
        if batch_key != "batch":
            adata_work.obs[batch_key] = adata_work.obs["batch"]
        SCVI.prepare_query_anndata(adata_work, reference_model_dir)
        model = SCVI.load_query_data(adata_work, reference_model_dir)
        plan_kwargs["weight_decay"] = 0.0
    else:
        raise ValueError(f"Invalid reference model type: {reference_model_type}")
else:
    categorical_covariates = "${categorical_covariates}"
    continuous_covariates = "${continuous_covariates}"

    categorical_covariates = categorical_covariates.split(",") if categorical_covariates else None
    continuous_covariates = continuous_covariates.split(",") if continuous_covariates else None

    SCVI.setup_anndata(
        adata_work,
        batch_key="${batch_col}",
        categorical_covariate_keys=categorical_covariates,
        continuous_covariate_keys=continuous_covariates,
    )

    model = SCVI(
        adata_work,
        n_hidden=int("${n_hidden}"),
        n_layers=int("${n_layers}"),
        n_latent=int("${n_latent}"),
        dispersion="${dispersion}",
        gene_likelihood="${gene_likelihood}",
        use_observed_lib_size="${use_observed_lib_size}" == "true",
    )

if "${task.accelerator ? 'true' : 'false'}" == "true":
    model.to_device(0)

model.train(
    early_stopping=True,
    max_epochs=int("${max_epochs}") if "${max_epochs?:''}" else None,
    plan_kwargs=plan_kwargs,
)

# Round to ensure hashes are stable
adata.obsm["X_emb"] = model.get_latent_representation()

write_anndata(adata, "${prefix}", "${anndata_format}")
model.save("${prefix}_model")

df = pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names)
df.to_pickle("X_${prefix}.pkl")

# Versions

versions = {"${task.process}": {"scvi": scvi.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
